#ifdef POLYSOLVE_WITH_TRILINOS

////////////////////////////////////////////////////////////////////////////////
#include <polysolve/LinearSolverTrilinos.hpp>
#include <string>
////////////////////////////////////////////////////////////////////////////////

namespace polysolve
{
    LinearSolverTrilinos::LinearSolverTrilinos()
    {
        precond_num_ = 0;
#ifdef HAVE_MPI
        int done_already;
        MPI_Initialized(&done_already);
        if (!done_already)
        {
            /* Initialize MPI */
            int argc = 1;
            char name[] = "";
            char *argv[] = {name};
            char **argvv = &argv[0];
            MPI_Init(&argc, &argvv);
            CommPtr = new Epetra_MpiComm(MPI_COMM_WORLD);
        }
#else
     CommPtr=new Epetra_SerialComm;
#endif
    }

    ////////////////////////////////////////////////////////////////
    void LinearSolverTrilinos::setParameters(const json &params)
    {
        if (params.contains("Trilinos"))
        {
            if (params["Trilinos"].contains("block_size"))
            {
                if (params["Trilinos"]["block_size"] == 2 || params["Trilinos"]["block_size"] == 3)
                {
                    numPDEs = params["Trilinos"]["block_size"];
                }
            }
            if (params["Trilinos"].contains("max_iter"))
            {
                max_iter_ = params["Trilinos"]["max_iter"];
            }
            if (params["Trilinos"].contains("tolerance"))
            {
                conv_tol_ = params["Trilinos"]["tolerance"];
            }
        }
    }

    /////////////////////////////////////////////////
    void LinearSolverTrilinos::getInfo(json &params) const
    {
        params["num_iterations"] = iterations_;
        params["final_res_norm"] = residual_error_;
    }

    /////////////////////////////////////////////////
    void LinearSolverTrilinos::factorize(const StiffnessMatrix &Ain)
    {
        assert(precond_num_ > 0);

        Eigen::SparseMatrix<double,Eigen::RowMajor> Arow(Ain);
        int mypid = CommPtr->MyPID();        
        int indexBase=0;
        int numGlobalElements = Arow.nonZeros();
        int numGlobalRows=Arow.rows();
        Epetra_Map *rowMap=NULL;
        int numNodes= numGlobalRows /numPDEs;       
        if ((numGlobalRows - numNodes * numPDEs) != 0 && !mypid){
            throw std::runtime_error("Number of matrix rows is not divisible by #dofs");
        }
        int numMyNodes;
        int nproc = CommPtr->NumProc();
        if (CommPtr->MyPID() < nproc-1) numMyNodes = numNodes / nproc;
        else numMyNodes = numNodes - (numNodes/nproc) * (nproc-1);
        rowMap = new Epetra_Map(numGlobalRows,numMyNodes*numPDEs,indexBase,(*CommPtr));

        A = new Epetra_CrsMatrix(Copy,*rowMap,0); //Can allocate memory for each row in advance

        {
            int nnzs=0;
            for (int k=0 ; k < Arow.outerSize(); k++)
            {
                // std::cout<<Arow.innerVector(k).nonZeros()<<" ";
                int numEntries=Arow.outerIndexPtr()[k+1]-Arow.outerIndexPtr()[k];
                int* indices=Arow.innerIndexPtr () + nnzs;
                double* values=Arow.valuePtr () + nnzs;
                A->InsertGlobalValues(k,numEntries,values,indices);
                nnzs=nnzs+numEntries;
            }
        }
        A->FillComplete();
    }

    namespace
    {
        void TrilinosML_SetDefaultOptions(Teuchos::ParameterList &MLList)
        {
            std::string aggr_type="MIS";
            double str_connect=0.08;
            ML_Epetra::SetDefaults("SA", MLList);
            MLList.set("aggregation: type",aggr_type); // Aggresive Method  
            // MLList.set("aggregation: type","Uncoupled"); // Fixed size
            

            MLList.set("aggregation: threshold",str_connect);
            
            // Smoother Settings
            MLList.set("smoother: type","Chebyshev");
            MLList.set("smoother: sweeps", 5); //Chebyshev degree
            MLList.set("smoother: Chbeshev alpha",30.0);

            //Coarser Settings
            MLList.set("coarse: max size",5);

            MLList.set("ML output",10);
        }
    }


    void LinearSolverTrilinos::solve(const Eigen::Ref<const VectorXd> rhs, Eigen::Ref<VectorXd> result)
    {
        int output=10; //how often to print residual history
        Teuchos::ParameterList MLList;
        TrilinosML_SetDefaultOptions(MLList);
        MLList.set("PDE equations",numPDEs);  
        ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

        Epetra_Vector x(A->RowMap());
        Epetra_Vector b(A->RowMap());
        for (size_t i = 0; i < rhs.size(); i++)
        {
            x[i]=result[i];
            b[i]=rhs[i];
        }
        // std::cout<<"x[0] "<<x[0]<<std::endl;
        // std::cout<<"b[0] "<<b[0]<<std::endl;

        Epetra_LinearProblem Problem(A,&x,&b);
        AztecOO solver(Problem);
        solver.SetAztecOption(AZ_solver, AZ_cg);
        solver.SetPrecOperator(MLPrec);
        solver.SetAztecOption(AZ_output, output);

        solver.Iterate(max_iter_, conv_tol_ );

        //Calculate a final residual
        Epetra_Vector workvec(A->RowMap());
        double mynorm;
        A->Multiply(false,x,workvec);
        workvec.Update(1.0,b,-1.0);
        b.Norm2(&mynorm);
        workvec.Scale(1./mynorm);
        workvec.Norm2(&mynorm);
        if (CommPtr->MyPID() == 0)
        {
            
            // std::cout<<"Max iterations "<<max_iter_<<std::endl;
            // std::cout<<"Residual is "<<mynorm<<std::endl;
            // std::cout<<"Trilinos Residual is "<<solver.ScaledResidual ()<<std::endl;
            // std::cout<<"Iterations are "<<solver.NumIters()<<std::endl;            
            residual_error_=mynorm;
            iterations_=solver.NumIters();
        }

        for (size_t i = 0; i < rhs.size(); i++)
        {
            result[i]=x[i];
        }
        delete MLPrec;
    }
    LinearSolverTrilinos:: ~LinearSolverTrilinos()
    {
        
        delete A;
#ifdef HAVE_MPI
        MPI_Finalize() ;
#endif
    }
}
  
#endif
