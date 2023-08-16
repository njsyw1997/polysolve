//////////////////////////////////////////////////////////////////////////
#include <polysolve/FEMSolver.hpp>
#include <catch2/catch.hpp>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <polysolve/LinearSolverAMGCL.hpp>
#include <polysolve/LinearSolverEigen.hpp>
#include <Eigen/PardisoSupport>
#include <polysolve/LinearSolver.hpp>
#include <polysolve/LinearSolverTrilinos.hpp>
#include <amgcl/io/mm.hpp>
//////////////////////////////////////////////////////////////////////////

using namespace polysolve;

void loadSymmetric(Eigen::SparseMatrix<double> &A, std::string PATH)
{
    std::ifstream fin(PATH);
    long int M, N, L;
    while (fin.peek() == '%')
    {
        fin.ignore(2048, '\n');
    }
    fin >> M >> N >> L;
    A.resize(M, N);
    A.reserve(L * 2 - M);
    std::vector<Eigen::Triplet<double>> triple;
    for (size_t i = 0; i < L; i++)
    {
        int m, n;
        double data;
        fin >> m >> n >> data;
        triple.push_back(Eigen::Triplet<double>(m - 1, n - 1, data));
        if (m != n)
        {
            triple.push_back(Eigen::Triplet<double>(n - 1, m - 1, data));
        }
    }
    fin.close();
    A.setFromTriplets(triple.begin(), triple.end());
};

void loadVec(Eigen::VectorXd &A, std::string PATH)
{
    std::ifstream fin(PATH);
    long int M, N;
    while (fin.peek() == '%')
    {
        fin.ignore(2048, '\n');
    }
    fin >> M >> N;
    A.resize(M*N);
    for (size_t i = 0; i < M*N; i++)
    {
        double data;
        fin >> data;
        A[i]=data;
    }
    fin.close();
};

// #ifdef POLYSOLVE_WITH_TRILINOS
// TEST_CASE("trilinos", "[solver]")
// {
//     std::cout << "Polysolve AMGCL Solver" << std::endl;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/stiffness.mtx");
//     REQUIRE(ok);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     loadVec(b, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/rhs.mtx");
//     Eigen::VectorXd points;
//     loadVec(points, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/vec.mtx");
//     test_vertices.resize(points.size()/3,3);
//     for (int i = 0; i < points.size()/3; i++)
//         for (int j = 0; j < 3; j++)
//             {
//                 test_vertices(i,j)=points(i+j*points.size()/3);
//             }
//     Eigen::VectorXd x_b(A.rows());
//     x_b.setZero();
//     Eigen::VectorXd x(A.rows());
//     x.setZero();
//     {
//         json solver_info;
//         auto solver = new LinearSolverTrilinos();

//         json params=R"({
//             "Trilinos": {
//                 "max_iter": 1000,
//                 "tolerance": 1e-8,
//                 "block_size": 3
//             }
//             })"_json;
//         params["Trilinos"]["block_size"] = 3;
//         // params["AMGCL"]["solver_type"] = "cg";
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         solver->getInfo(solver_info);
//         REQUIRE(solver_info["num_iterations"] > 0);
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//         delete solver;
//     }
//     REQUIRE((A * x_b - b).norm() / b.norm() < 1e-7);
// }
// #endif

#ifdef POLYSOLVE_WITH_TRILINOS
TEST_CASE("trilinos", "[solver]")
{
    std::cout << "Polysolve AMGCL Solver" << std::endl;
    Eigen::SparseMatrix<double> A;
    const bool ok = loadMarket(A, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/stiffness.mtx");
    REQUIRE(ok);
    std::cout << "Matrix Load OK" << std::endl;
    Eigen::VectorXd b(A.rows());
    loadVec(b, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/rhs.mtx");
    Eigen::VectorXd points;
    loadVec(points, "/home/yiwei/test_bar/tmp/Trilinos/square_beam_1.25/bar_sa/P1/ref1/block3/Thread8/0/output/vec.mtx");
    test_vertices.resize(points.size()/3,3);
    for (int i = 0; i < points.size()/3; i++)
        for (int j = 0; j < 3; j++)
            {
                test_vertices(i,j)=points(i+j*points.size()/3);
            }
    Eigen::VectorXd x_b(A.rows());
    x_b.setZero();
    Eigen::VectorXd x(A.rows());
    x.setZero();
    {
        json solver_info;
        auto solver = new LinearSolverEigenDirect<Eigen::PardisoLLT<polysolve::StiffnessMatrix>>("Eigen::PardisoLLT");

        json params=R"({
            "Pardiso": {
                "mtype": 1
            }
            })"_json;
        // params["AMGCL"]["solver_type"] = "cg";
        solver->setParameters(params);
        solver->analyzePattern(A, A.rows());
        solver->factorize(A);
        solver->solve(b, x_b);
        solver->getInfo(solver_info);
        // REQUIRE(solver_info["num_iterations"] > 0);
        std::cout << solver_info["num_iterations"] << std::endl;
        std::cout << solver_info["final_res_norm"] << std::endl;
        delete solver;
    }
    REQUIRE((A * x_b - b).norm() / b.norm() < 1e-7);
}
#endif

