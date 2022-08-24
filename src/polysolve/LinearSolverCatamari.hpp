#pragma once

#ifdef POLYSOLVE_WITH_CATAMARI

////////////////////////////////////////////////////////////////////////////////
#include <polysolve/LinearSolver.hpp>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include "catamari.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// https://computation.llnl.gov/sites/default/files/public/hypre-2.11.2_usr_manual.pdf
// https://github.com/LLNL/hypre/blob/v2.14.0/docs/HYPRE_usr_manual.pdf
//

namespace polysolve
{

    class LinearSolverCatamari : public LinearSolver
    {

    public:
        LinearSolverCatamari();
        ~LinearSolverCatamari();

    private:
        POLYSOLVE_DELETE_MOVE_COPY(LinearSolverCatamari)

    public:
        //////////////////////
        // Public interface //
        //////////////////////

        // Set solver parameters
        virtual void setParameters(const json &params) override;

        // Retrieve memory information from Pardiso
        virtual void getInfo(json &params) const override;

        // Analyze sparsity pattern
        virtual void analyzePattern(const StiffnessMatrix &A, const int precond_num) override { precond_num_ = precond_num; }

        // Factorize system matrix
        virtual void factorize(const StiffnessMatrix &A) override;

        // Solve the linear system Ax = b
        virtual void solve(const Ref<const VectorXd> b, Ref<VectorXd> x) override;

        // Name of the solver type (for debugging purposes)
        virtual std::string name() const override { return "Catamari"; }

    private:
        int precond_num_;
        catamari::CoordinateMatrix<double> A;
        catamari::SparseLDL<double> ldl;
        catamari::SparseLDLControl<double> ldl_control;
        catamari::SparseLDLResult<double> factor;
        catamari::BlasMatrix<double> right_hand_sides;
    };

} // namespace polysolve

#endif
