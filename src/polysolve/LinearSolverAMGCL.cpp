#ifdef POLYSOLVE_WITH_AMGCL

////////////////////////////////////////////////////////////////////////////////
#include <polysolve/LinearSolverAMGCL.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace polysolve
{
    bool is_test = true;
    namespace
    {
        /* https://stackoverflow.com/questions/15904896/range-based-for-loop-on-a-dynamic-array */
        template <typename T>
        struct WrappedArray
        {
            WrappedArray(const T *first, const T *last)
                : begin_{first}, end_{last} {}
            WrappedArray(const T *first, std::ptrdiff_t size)
                : WrappedArray{first, first + size} {}

            const T *begin() const noexcept { return begin_; }
            const T *end() const noexcept { return end_; }
            const T &operator[](const size_t i) const { return begin_[i]; }

            const T *begin_;
            const T *end_;
        };

        json default_params()
        {
            json params = R"({
            })"_json;

            return params;
        }

        void set_params(const json &params, json &out)
        {
            if (params.contains("AMGCL"))
            {
                // Patch the stored params with input ones
                // if (params["AMGCL"].contains("precond"))
                //     out["precond"].merge_patch(params["AMGCL"]["precond"]);
                // if (params["AMGCL"].contains("solver"))
                //     out["solver"].merge_patch(params["AMGCL"]["solver"]);
                out = params["AMGCL"];
                // if (out["precond"]["class"] == "schur_pressure_correction")
                // {
                //     // Initialize the u and p solvers with a tolerance that is comparable to the main solver's
                //     if (!out["precond"].contains("usolver"))
                //     {
                //         out["precond"]["usolver"] = R"({"solver": {"maxiter": 100}})"_json;
                //         out["precond"]["usolver"]["solver"]["tol"] = 10 * out["solver"]["tol"].get<double>();
                //     }
                //     if (!out["precond"].contains("usolver"))
                //     {
                //         out["precond"]["psolver"] = R"({"solver": {"maxiter": 100}})"_json;
                //         out["precond"]["psolver"]["solver"]["tol"] = 10 * out["solver"]["tol"].get<double>();
                //     }
                // }
            }
        }
    } // namespace

    ////////////////////////////////////////////////////////////////////////////////

    LinearSolverAMGCL::LinearSolverAMGCL()
    {
        // NOTE: usolver and psolver parameters are only used if the
        // preconditioner class is "schur_pressure_correction"
        precond_num_ = 0;
    }

    // Set solver parameters
    void LinearSolverAMGCL::setParameters(const json &params)
    {
        if (params.contains("AMGCL"))
        {
            // Specially named parameters to match other solvers
            if (params["AMGCL"].contains("block_size"))
            {
                block_size_ = params["AMGCL"]["block_size"];
            }
            if (block_size_ == 2)
            {
                block2_solver_.setParameters(params);
                return;
            }
            else if (block_size_ == 3)
            {
                block3_solver_.setParameters(params);
                return;
            }

            set_params(params, params_);
        }
    }

    void LinearSolverAMGCL::getInfo(json &params) const
    {
        if (block_size_ == 2)
        {
            block2_solver_.getInfo(params);
            return;
        }
        else if (block_size_ == 3)
        {
            block3_solver_.getInfo(params);
            return;
        }
        params["num_iterations"] = iterations_;
        params["final_res_norm"] = residual_error_;
    }

    ////////////////////////////////////////////////////////////////////////////////

    void LinearSolverAMGCL::factorize(const StiffnessMatrix &Ain)
    {
        if (block_size_ == 2)
        {
            block2_solver_.factorize(Ain);
            return;
        }
        else if (block_size_ == 3)
        {
            block3_solver_.factorize(Ain);
            return;
        }
        assert(precond_num_ > 0);

        int numRows = Ain.rows();

        WrappedArray<StiffnessMatrix::StorageIndex> ia(Ain.outerIndexPtr(), numRows + 1);
        WrappedArray<StiffnessMatrix::StorageIndex> ja(Ain.innerIndexPtr(), Ain.nonZeros());
        WrappedArray<StiffnessMatrix::Scalar> a(Ain.valuePtr(), Ain.nonZeros());
        if (params_["precond"]["class"] == "schur_pressure_correction")
        {
            std::vector<char> pmask(numRows, 0);
            for (size_t i = precond_num_; i < numRows; ++i)
                pmask[i] = 1;
            params_["precond"]["pmask"] = pmask;
        }

        // AMGCL takes the parameters as a Boost property_tree (i.e., another JSON data structure)
        std::stringstream ss_params;
        ss_params << params_;
        boost::property_tree::ptree pt_params;
        boost::property_tree::read_json(ss_params, pt_params);
        auto A = std::tie(numRows, ia, ja, a);
        prof.tic("setup");
        solver_ = std::make_unique<Solver>(A, pt_params);
        prof.toc("setup");
        iterations_ = 0;
        residual_error_ = 0;
    }

    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////

    void LinearSolverAMGCL::solve(const Eigen::Ref<const VectorXd> rhs, Eigen::Ref<VectorXd> result)
    {
        if (block_size_ == 2)
        {
            block2_solver_.solve(rhs, result);
            return;
        }
        else if (block_size_ == 3)
        {
            block3_solver_.solve(rhs, result);
            return;
        }
        assert(result.size() == rhs.size());
        std::vector<double> _rhs(rhs.data(), rhs.data() + rhs.size());
        if (is_test)
        {
            result=Eigen::VectorXd::Random(result.size());
        }

        std::vector<double> x(result.data(), result.data() + result.size());
        auto rhs_b = Backend::copy_vector(_rhs, backend_params_);
        auto x_b = Backend::copy_vector(x, backend_params_);

        assert(solver_ != nullptr);
        prof.tic("solve");
        std::tie(iterations_, residual_error_) = (*solver_)(*rhs_b, *x_b);
        prof.toc("solve");
            std::copy(&(*x_b)[0], &(*x_b)[0] + result.size(), result.data());
        if (is_test)
        {
            std::cout << (*solver_) << std::endl;
            std::cout << prof << std::endl;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    LinearSolverAMGCL::~LinearSolverAMGCL()
    {
    }

    template <int BLOCK_SIZE>
    LinearSolverAMGCL_Block<BLOCK_SIZE>::LinearSolverAMGCL_Block()
    {
        // NOTE: usolver and psolver parameters are only used if the
        // preconditioner class is "schur_pressure_correction"
        precond_num_ = 0;
    }

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_Block<BLOCK_SIZE>::setParameters(const json &params)
    {
        set_params(params, params_);
    }

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_Block<BLOCK_SIZE>::getInfo(json &params) const
    {
        params["num_iterations"] = iterations_;
        params["final_res_norm"] = residual_error_;
    }

    ////////////////////////////////////////////////////////////////////////////////

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_Block<BLOCK_SIZE>::factorize(const StiffnessMatrix &Ain)
    {
        assert(precond_num_ > 0);

        int numRows = Ain.rows();

        WrappedArray<StiffnessMatrix::StorageIndex> ia(Ain.outerIndexPtr(), numRows + 1);
        WrappedArray<StiffnessMatrix::StorageIndex> ja(Ain.innerIndexPtr(), Ain.nonZeros());
        WrappedArray<StiffnessMatrix::Scalar> a(Ain.valuePtr(), Ain.nonZeros());

        if (params_["precond"]["class"] == "schur_pressure_correction")
        {
            std::vector<char> pmask(numRows, 0);
            for (size_t i = precond_num_; i < numRows; ++i)
                pmask[i] = 1;
            params_["precond"]["pmask"] = pmask;
        }

        // AMGCL takes the parameters as a Boost property_tree (i.e., another JSON data structure)
        std::stringstream ss_params;
        ss_params << params_;
        boost::property_tree::ptree pt_params;
        boost::property_tree::read_json(ss_params, pt_params);

        auto A = std::tie(numRows, ia, ja, a);
        auto Ab = amgcl::adapter::block_matrix<dmat_type>(A);
        prof.tic("setup");
        solver_ = std::make_unique<Solver>(Ab, pt_params);
        prof.toc("setup");
        iterations_ = 0;
        residual_error_ = 0;
    }

    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_Block<BLOCK_SIZE>::solve(const Eigen::Ref<const VectorXd> rhs, Eigen::Ref<VectorXd> result)
    {
        assert(result.size() == rhs.size());
        std::vector<double> _rhs(rhs.data(), rhs.data() + rhs.size());
        // if (is_test)
        // {
        //     result=Eigen::VectorXd::Random(result.size());
        // }
        std::vector<double> x(result.data(), result.data() + result.size());
        auto rhs_b = amgcl::backend::reinterpret_as_rhs<dmat_type>(_rhs);
        auto x_b = amgcl::backend::reinterpret_as_rhs<dmat_type>(x);

        assert(solver_ != nullptr);
        prof.tic("solve");
        std::tie(iterations_, residual_error_) = (*solver_)(rhs_b, x_b);
        prof.toc("solve");
        for (size_t i = 0; i < rhs.size() / BLOCK_SIZE; i++)
            for (size_t j = 0; j < BLOCK_SIZE; j++)
            {
                result[BLOCK_SIZE * i + j] = x_b[i](j);
            }
        if (is_test)
        {
            std::cout << (*solver_) << std::endl;
            std::cout << prof << std::endl;
        }
    }

    template <int BLOCK_SIZE>
    LinearSolverAMGCL_Block<BLOCK_SIZE>::~LinearSolverAMGCL_Block()
    {
    }


    template <int BLOCK_SIZE>
    LinearSolverAMGCL_NS<BLOCK_SIZE>::LinearSolverAMGCL_NS()
    {
        // NOTE: usolver and psolver parameters are only used if the
        // preconditioner class is "schur_pressure_correction"
        precond_num_ = 0;
    }

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_NS<BLOCK_SIZE>::setParameters(const json &params)
    {
        set_params(params, params_);
    }

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_NS<BLOCK_SIZE>::getInfo(json &params) const
    {
        params["num_iterations"] = iterations_;
        params["final_res_norm"] = residual_error_;
    }

    ////////////////////////////////////////////////////////////////////////////////

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_NS<BLOCK_SIZE>::factorize(const StiffnessMatrix &Ain)
    {
        assert(precond_num_ > 0);

        int numRows = Ain.rows();

        WrappedArray<StiffnessMatrix::StorageIndex> ia(Ain.outerIndexPtr(), numRows + 1);
        WrappedArray<StiffnessMatrix::StorageIndex> ja(Ain.innerIndexPtr(), Ain.nonZeros());
        WrappedArray<StiffnessMatrix::Scalar> a(Ain.valuePtr(), Ain.nonZeros());

        if (params_["precond"]["class"] == "schur_pressure_correction")
        {
            std::vector<char> pmask(numRows, 0);
            for (size_t i = precond_num_; i < numRows; ++i)
                pmask[i] = 1;
            params_["precond"]["pmask"] = pmask;
        }

        // AMGCL takes the parameters as a Boost property_tree (i.e., another JSON data structure)
        std::stringstream ss_params;
        ss_params << params_;
        boost::property_tree::ptree pt_params;
        boost::property_tree::read_json(ss_params, pt_params);

        auto A = std::tie(numRows, ia, ja, a);
        auto Ab = amgcl::adapter::block_matrix<dmat_type>(A);
        prof.tic("setup");
        solver_ = std::make_unique<Solver>(Ab, pt_params);
        prof.toc("setup");
        iterations_ = 0;
        residual_error_ = 0;
    }

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_NS<BLOCK_SIZE>::factorize(const StiffnessMatrix &Ain, const VectorXd &coo)
    {
        assert(precond_num_ > 0);

        int numRows = Ain.rows();

        WrappedArray<StiffnessMatrix::StorageIndex> ia(Ain.outerIndexPtr(), numRows + 1);
        WrappedArray<StiffnessMatrix::StorageIndex> ja(Ain.innerIndexPtr(), Ain.nonZeros());
        WrappedArray<StiffnessMatrix::Scalar> a(Ain.valuePtr(), Ain.nonZeros());

        if (params_["precond"]["class"] == "schur_pressure_correction")
        {
            std::vector<char> pmask(numRows, 0);
            for (size_t i = precond_num_; i < numRows; ++i)
                pmask[i] = 1;
            params_["precond"]["pmask"] = pmask;
        }

        // AMGCL takes the parameters as a Boost property_tree (i.e., another JSON data structure)
        std::stringstream ss_params;
        ss_params << params_;
        boost::property_tree::ptree pt_params;
        boost::property_tree::read_json(ss_params, pt_params);
        typename Solver::params prm(pt_params);

        prm.precond.coarsening.nullspace.cols = amgcl::coarsening::rigid_body_modes(
        BLOCK_SIZE, coo, prm.precond.coarsening.nullspace.B);

        auto A = std::tie(numRows, ia, ja, a);
        auto Ab = amgcl::adapter::block_matrix<dmat_type>(A);
        prof.tic("setup");
        solver_ = std::make_unique<Solver>(Ab, prm);
        prof.toc("setup");
        iterations_ = 0;
        residual_error_ = 0;
    }

    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////

    template <int BLOCK_SIZE>
    void LinearSolverAMGCL_NS<BLOCK_SIZE>::solve(const Eigen::Ref<const VectorXd> rhs, Eigen::Ref<VectorXd> result)
    {
        assert(result.size() == rhs.size());
        std::vector<double> _rhs(rhs.data(), rhs.data() + rhs.size());
        // if (is_test)
        // {
        //     result=Eigen::VectorXd::Random(result.size());
        // }
        std::vector<double> x(result.data(), result.data() + result.size());
        auto rhs_b = amgcl::backend::reinterpret_as_rhs<dmat_type>(_rhs);
        auto x_b = amgcl::backend::reinterpret_as_rhs<dmat_type>(x);

        assert(solver_ != nullptr);
        prof.tic("solve");
        std::tie(iterations_, residual_error_) = (*solver_)(rhs_b, x_b);
        prof.toc("solve");
        for (size_t i = 0; i < rhs.size() / BLOCK_SIZE; i++)
            for (size_t j = 0; j < BLOCK_SIZE; j++)
            {
                result[BLOCK_SIZE * i + j] = x_b[i](j);
            }
        if (is_test)
        {
            std::cout << (*solver_) << std::endl;
            std::cout << prof << std::endl;
        }
    }

    template <int BLOCK_SIZE>
    LinearSolverAMGCL_NS<BLOCK_SIZE>::~LinearSolverAMGCL_NS()
    {
    }


    template class LinearSolverAMGCL_Block<2>;
    template class LinearSolverAMGCL_Block<3>;
    template class LinearSolverAMGCL_NS<2>;
    template class LinearSolverAMGCL_NS<3>;
} // namespace polysolve

#endif
