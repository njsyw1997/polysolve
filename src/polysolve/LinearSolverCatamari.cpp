#ifdef POLYSOLVE_WITH_CATAMARI

////////////////////////////////////////////////////////////////////////////////
#include <polysolve/LinearSolverCatamari.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
////////////////////////////////////////////////////////////////////////////////
namespace polysolve
{
    LinearSolverCatamari::LinearSolverCatamari()
    {
        precond_num_ = 0;
    }

    void LinearSolverCatamari::setParameters(const json &params)
    {
        /*
        Options for factorization:
        catamari::kCholeskyFactorization,
        catamari::kLDLAdjointFactorization,
        catamari::kLDLTransposeFactorization.
        */
        ldl_control.SetFactorizationType(catamari::kCholeskyFactorization);
    }

    void LinearSolverCatamari::getInfo(json &params) const
    {
        // Catamari returns nothing
    }

    void LinearSolverCatamari::factorize(const StiffnessMatrix &Ain)
    {
        // Set up matrix
        A.Resize(Ain.rows(), Ain.cols());
        A.ReserveEntryAdditions(Ain.nonZeros());
        for (int k = 0; k < Ain.outerSize(); ++k)
        {
            for (StiffnessMatrix::InnerIterator it(Ain, k); it; ++it)
            {
                A.QueueEntryAddition(it.row(), it.col(), it.value());
            }
        }
        A.FlushEntryQueues();

        // Factorize the matrix
        factor = ldl.Factor(A, ldl_control);
    }

    void LinearSolverCatamari::solve(const Eigen::Ref<const VectorXd> rhs, Eigen::Ref<VectorXd> result)
    {
        // Copy rhs
        right_hand_sides.Resize(rhs.rows(), rhs.cols());
        // The (i, j) entry of the right-hand side can easily be read or modified, e.g.:
        for (int i = 0; i < rhs.rows(); i++)
        {
            for (int j = 0; j < rhs.cols(); j++)
            {
                right_hand_sides(i, j) = rhs(i, j);
            }
        }
        // Solve a linear system using the factorization.
        ldl.Solve(&right_hand_sides.view);
        // copy solution
        result.resize(rhs.rows(), rhs.cols());
        for (int i = 0; i < rhs.rows(); i++)
        {
            for (int j = 0; j < rhs.cols(); j++)
            {
                result(i, j) = right_hand_sides(i, j);
            }
        }
    }
    LinearSolverCatamari::~LinearSolverCatamari(){}

} // namespace polysolve

#endif