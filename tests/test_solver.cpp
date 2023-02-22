//////////////////////////////////////////////////////////////////////////
#include <polysolve/FEMSolver.hpp>
#include <catch2/catch.hpp>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <fstream>
#include <vector>
#include <ctime>
#include <polysolve/LinearSolverAMGCL.hpp>
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
    A.resize(M);
    for (size_t i = 0; i < M; i++)
    {
        double data;
        fin >> data;
        A[i]=data;
    }
    fin.close();
};
// TEST_CASE("all", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     auto solvers = LinearSolver::availableSolvers();

//     for (const auto &s : solvers)
//     {
//         if (s == "Eigen::DGMRES")
//             continue;
// #ifdef WIN32
//         if (s == "Eigen::ConjugateGradient" || s == "Eigen::BiCGSTAB" || s == "Eigen::GMRES" || s == "Eigen::MINRES")
//             continue;
// #endif
//         auto solver = LinearSolver::create(s, "");
//         json params;
//         params[s]["tolerance"] = 1e-10;
//         solver->setParameters(params);
//         Eigen::VectorXd b(A.rows());
//         b.setRandom();
//         Eigen::VectorXd x(b.size());
//         x.setZero();

//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);

//         // solver->getInfo(solver_info);

//         // std::cout<<"Solver error: "<<x<<std::endl;
//         const double err = (A * x - b).norm();
//         INFO("solver: " + s);
//         REQUIRE(err < 1e-8);
//     }
// }

// TEST_CASE("eigen_params", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     auto solvers = LinearSolver::availableSolvers();

//     for (const auto &s : solvers)
//     {
//         if (s == "Eigen::ConjugateGradient" || s == "Eigen::BiCGSTAB" || s == "Eigen::GMRES" || s == "Eigen::MINRES" || s == "Eigen::LeastSquaresConjugateGradient" || s == "Eigen::DGMRES")
//         {
//             auto solver = LinearSolver::create(s, "");
//             json params;
//             params[s]["max_iter"] = 1000;
//             params[s]["tolerance"] = 1e-10;
//             solver->setParameters(params);

//             Eigen::VectorXd b(A.rows());
//             b.setRandom();
//             Eigen::VectorXd x(b.size());
//             x.setZero();

//             solver->analyzePattern(A, A.rows());
//             solver->factorize(A);
//             solver->solve(b, x);

//             // solver->getInfo(solver_info);

//             // std::cout<<"Solver error: "<<x<<std::endl;
//             const double err = (A * x - b).norm();
//             INFO("solver: " + s);
//             REQUIRE(err < 1e-8);
//         }
//     }
// }

// TEST_CASE("pre_factor", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     auto solvers = LinearSolver::availableSolvers();

//     for (const auto &s : solvers)
//     {
//         if (s == "Eigen::DGMRES")
//             continue;
// #ifdef WIN32
//         if (s == "Eigen::ConjugateGradient" || s == "Eigen::BiCGSTAB" || s == "Eigen::GMRES" || s == "Eigen::MINRES")
//             continue;
// #endif
//         auto solver = LinearSolver::create(s, "");
//         solver->analyzePattern(A, A.rows());

//         std::default_random_engine eng{42};
//         std::uniform_real_distribution<double> urd(0.1, 5);

//         for (int i = 0; i < 10; ++i)
//         {
//             std::vector<Eigen::Triplet<double>> tripletList;

//             for (int k = 0; k < A.outerSize(); ++k)
//             {
//                 for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
//                 {
//                     if (it.row() == it.col())
//                     {
//                         tripletList.emplace_back(it.row(), it.col(), urd(eng) * 100);
//                     }
//                     else if (it.row() < it.col())
//                     {
//                         const double val = -urd(eng);
//                         tripletList.emplace_back(it.row(), it.col(), val);
//                         tripletList.emplace_back(it.col(), it.row(), val);
//                     }
//                 }
//             }

//             Eigen::SparseMatrix<double> Atmp(A.rows(), A.cols());
//             Atmp.setFromTriplets(tripletList.begin(), tripletList.end());

//             Eigen::VectorXd b(Atmp.rows());
//             b.setRandom();
//             Eigen::VectorXd x(b.size());
//             x.setZero();

//             solver->factorize(Atmp);
//             solver->solve(b, x);

//             // solver->getInfo(solver_info);

//             // std::cout<<"Solver error: "<<x<<std::endl;
//             const double err = (Atmp * x - b).norm();
//             INFO("solver: " + s);
//             REQUIRE(err < 1e-8);
//         }
//     }
// }

// #ifdef POLYSOLVE_WITH_HYPRE
// TEST_CASE("hypre", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     auto solver = LinearSolver::create("Hypre", "");
//     // solver->setParameters(params);
//     Eigen::VectorXd b(A.rows());
//     b.setRandom();
//     Eigen::VectorXd x(b.size());
//     x.setZero();

//     solver->analyzePattern(A, A.rows());
//     solver->factorize(A);
//     solver->solve(b, x);

//     // solver->getInfo(solver_info);

//     // std::cout<<"Solver error: "<<x<<std::endl;
//     const double err = (A * x - b).norm();
//     REQUIRE(err < 1e-8);
// }

// TEST_CASE("hypre_initial_guess", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     // solver->setParameters(params);
//     Eigen::VectorXd b(A.rows());
//     b.setRandom();
//     Eigen::VectorXd x(A.rows());
//     x.setZero();
//     {
//         json solver_info;

//         auto solver = LinearSolver::create("Hypre", "");
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] > 1);
//     }

//     {
//         json solver_info;
//         auto solver = LinearSolver::create("Hypre", "");
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);

//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] == 1);
//     }

//     // std::cout<<"Solver error: "<<x<<std::endl;
//     const double err = (A * x - b).norm();
//     REQUIRE(err < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_AMGCL
// TEST_CASE("amgcl_initial_guess", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     // solver->setParameters(params);
//     Eigen::VectorXd b(A.rows());
//     b.setRandom();
//     Eigen::VectorXd x(A.rows());
//     x.setZero();
//     {
//         json solver_info;

//         auto solver = LinearSolver::create("AMGCL", "");
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] > 0);
//     }

//     {
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);

//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] == 0);
//     }

//     // std::cout<<"Solver error: "<<x<<std::endl;
//     const double err = (A * x - b).norm();
//     REQUIRE(err < 1e-8);
// }
// #endif

// TEST_CASE("saddle_point_test", "[solver]")
// {
// #ifdef WIN32
// #ifndef NDEBUG
//     return;
// #endif
// #endif
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     bool ok = loadMarket(A, path + "/A0.mat");
//     REQUIRE(ok);

//     Eigen::VectorXd b;
//     ok = loadMarketVector(b, path + "/b0.mat");
//     REQUIRE(ok);

//     auto solver = LinearSolver::create("SaddlePointSolver", "");
//     solver->analyzePattern(A, 9934);
//     solver->factorize(A);
//     Eigen::VectorXd x(A.rows());
//     solver->solve(b, x);
//     const double err = (A * x - b).norm();
//     REQUIRE(err < 1e-8);
// }

// #ifdef POLYSOLVE_WITH_AMGCL
// TEST_CASE("amgcl_blocksolver_small_scale", "[solver]")
// {
// // #ifndef NDEBUG
// //     return;
// // #endif
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);

//     // solver->setParameters(params);
//     Eigen::VectorXd b(A.rows());
//     b.setRandom();
//     Eigen::VectorXd x(A.rows());
//     Eigen::VectorXd x_b(A.rows());
//     x.setZero();
//     x_b.setZero();
//     json params;
//     {
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] > 0);
//         const double err = (A * x - b).norm();
//         REQUIRE(err < 1e-5);
//     }

//     {
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         params["AMGCL"]["block_size"] = 3;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         solver->getInfo(solver_info);

//         REQUIRE(solver_info["num_iterations"] > 0);
//         const double err = (A * x_b - b).norm();
//         REQUIRE(err < 1e-5);
//     }
// }
// #endif

// #ifdef POLYSOLVE_WITH_AMGCL
// TEST_CASE("amgcl_blocksolver_b2", "[solver]")
// {
// #ifndef NDEBUG
//     return;
// #endif
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "gr_30_30.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);

//     std::cout << "Matrix Load OK" << std::endl;

//     Eigen::VectorXd b(A.rows());
//     b.setRandom();
//     Eigen::VectorXd x(A.rows());
//     Eigen::VectorXd x_b(A.rows());
//     x.setOnes();
//     x_b.setOnes();
//     {
//         amgcl::profiler<> prof("gr_30_30_Scalar");
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         prof.tic("setup");
//         json params;
//         params["AMGCL"]["tolerance"] = 1e-8;
//         params["AMGCL"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         prof.toc("setup");
//         prof.tic("solve");
//         solver->solve(b, x);
//         prof.toc("solve");
//         solver->getInfo(solver_info);
//         REQUIRE(solver_info["num_iterations"] > 0);
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl
//                   << prof << std::endl;
//     }
//     {
//         amgcl::profiler<> prof("gr_30_30_Block");
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         prof.tic("setup");
//         json params;
//         params["AMGCL"]["tolerance"] = 1e-8;
//         params["AMGCL"]["max_iter"] = 1000;
//         params["AMGCL"]["block_size"] = 2;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         prof.toc("setup");
//         prof.tic("solve");
//         solver->solve(b, x_b);
//         prof.toc("solve");
//         solver->getInfo(solver_info);
//         REQUIRE(solver_info["num_iterations"] > 0);
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl
//                   << prof << std::endl;
//     }
//     REQUIRE((A * x - b).norm() / b.norm() < 1e-7);
//     REQUIRE((A * x_b - b).norm() / b.norm() < 1e-7);
// }
// #endif

// #ifdef POLYSOLVE_WITH_AMGCL
// TEST_CASE("amgcl_blocksolver_crystm03_CG", "[solver]")
// {
// #ifndef NDEBUG
//     return;
// #endif
//     std::cout << "Polysolve AMGCL Solver" << std::endl;
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "crystm03.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x_b(A.rows());
//     x_b.setZero();
//     Eigen::VectorXd x(A.rows());
//     x.setZero();
//     {
//         amgcl::profiler<> prof("crystm03_Block");
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         prof.tic("setup");
//                 json params;
//         params["AMGCL"]["tolerance"] = 1e-8;
//         params["AMGCL"]["max_iter"] = 1000;
//         params["AMGCL"]["block_size"] = 3;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         prof.toc("setup");
//         prof.tic("solve");
//         solver->solve(b, x_b);
//         prof.toc("solve");
//         solver->getInfo(solver_info);
//         REQUIRE(solver_info["num_iterations"] > 0);
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl
//                   << prof << std::endl;
//     }
//     {
//         amgcl::profiler<> prof("crystm03_Scalar");
//         json solver_info;
//         auto solver = LinearSolver::create("AMGCL", "");
//         prof.tic("setup");
//                 json params;
//         params["AMGCL"]["tolerance"] = 1e-8;
//         params["AMGCL"]["max_iter"] = 10000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         prof.toc("setup");
//         prof.tic("solve");
//         solver->solve(b, x);
//         prof.toc("solve");
//         solver->getInfo(solver_info);
//         REQUIRE(solver_info["num_iterations"] > 0);
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl
//                   << prof << std::endl;
//     }
//     REQUIRE((A * x - b).norm() / b.norm() < 1e-7);
//     REQUIRE((A * x_b - b).norm() / b.norm() < 1e-7);
// }
// #endif

#ifdef POLYSOLVE_WITH_AMGCL
TEST_CASE("amgcl_blocksolver_crystm03_Bicgstab", "[solver]")
{
    std::cout << "Polysolve AMGCL Solver" << std::endl;
    Eigen::SparseMatrix<double> A;
    const bool ok = loadMarket(A, "/home/yiwei/result_newblock/AMGCL_0.08/square_beam_structed_25/bar_sa/P1/ref0/block3/Thread8/0/output/stiffness.mtx");
    REQUIRE(ok);
    std::cout << "Matrix Load OK" << std::endl;

    Eigen::VectorXd b(A.rows());
    loadVec(b, "/home/yiwei/result_newblock/AMGCL_0.08/square_beam_structed_25/bar_sa/P1/ref0/block3/Thread8/0/output/rhs.mtx");
    Eigen::VectorXd points;
    loadVec(points, "/home/yiwei/result_newblock/AMGCL_0.08/square_beam_structed_25/bar_sa/P1/ref0/block3/Thread8/0/output/points.mtx");
    Eigen::VectorXd x_b(A.rows());
    x_b.setZero();
    Eigen::VectorXd x(A.rows());
    x.setZero();
    {
        amgcl::profiler<> prof("crystm03_Block");
        json solver_info;
        auto solver = new LinearSolverAMGCL();
        prof.tic("setup");
        json params=R"({
        "AMGCL":{ 
        "precond": {
            "relax": {
                "degree": 5,
                "type": "chebyshev",
                "power_iters": 100,
                "higher": 1,
                "lower": 0.03,
                "scale": true
            },
            "class": "amg",
            "max_levels": 25,
            "direct_coarse": true,
            "ncycle": 1,
            "block_size": 3,
            "coarsening": {
                "type": "smoothed_aggregation",
                "estimate_spectral_radius": false,
                "relax": 1.0,
                "aggr": {
                    "eps_strong": 0.08
                }
            }
        },
        "solver": {
            "tol": 1e-8,
            "maxiter": 1000,
            "type": "cg",
            "verbose":true
        }
        }
        })"_json;
        params["AMGCL"]["tolerance"] = 1e-8;
        params["AMGCL"]["max_iter"] = 1000;
        params["AMGCL"]["block_size"] = 3;
        // params["AMGCL"]["solver_type"] = "cg";
        solver->setParameters(params);
        solver->analyzePattern(A, A.rows());
        solver->factorize(A);
        prof.toc("setup");
        prof.tic("solve");
        solver->solve(b, x_b);
        prof.toc("solve");
        solver->getInfo(solver_info);
        REQUIRE(solver_info["num_iterations"] > 0);
        std::cout << solver_info["num_iterations"] << std::endl;
        std::cout << solver_info["final_res_norm"] << std::endl
                  << prof << std::endl;
        delete solver;
    }
    {
        amgcl::profiler<> prof("crystm03_nullspace");
        json solver_info;
        auto solver = new LinearSolverAMGCL();
        prof.tic("setup");
        json params=R"({
        "AMGCL":{ 
        "precond": {
            "relax": {
                "degree": 5,
                "type": "chebyshev",
                "power_iters": 100,
                "higher": 1,
                "lower": 0.03,
                "scale": true
            },
            "class": "amg",
            "max_levels": 25,
            "direct_coarse": true,
            "ncycle": 1,
            "block_size": 3,
            "coarsening": {
                "type": "smoothed_aggregation",
                "estimate_spectral_radius": false,
                "relax": 1.0,
                "aggr": {
                    "eps_strong": 0.08
                }
            }
        },
        "solver": {
            "tol": 1e-8,
            "maxiter": 1000,
            "type": "cg",
            "verbose":true
        }
        }
        })"_json;
        params["AMGCL"]["tolerance"] = 1e-8;
        params["AMGCL"]["max_iter"] = 1000;
        params["AMGCL"]["block_size"] = 3;
        // params["AMGCL"]["solver_type"] = "bicgstab";
        solver->setParameters(params);
        solver->analyzePattern(A, A.rows());
        solver->factorize(A,points);
        prof.toc("setup");
        prof.tic("solve");
        solver->solve(b, x);
        prof.toc("solve");
        solver->getInfo(solver_info);
        REQUIRE(solver_info["num_iterations"] > 0);
        std::cout << solver_info["num_iterations"] << std::endl;
        std::cout << solver_info["final_res_norm"] << std::endl
                  << prof << std::endl;
        delete solver;
    }
    REQUIRE((A * x - b).norm() / b.norm() < 1e-7);
    REQUIRE((A * x_b - b).norm() / b.norm() < 1e-7);
}
#endif

// #ifdef POLYSOLVE_WITH_HYPRE
// TEST_CASE("Hyprel_b2", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "gr_30_30.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["block_size"] = 2;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_HYPRE
// TEST_CASE("Hybre_crystm03", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "crystm03.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["block_size"] = 3;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_HYPRE
// TEST_CASE("hypre_smallscale", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Hypre", "");
//         json params;
//         params["Hypre"]["block_size"] = 3;
//         params["Hypre"]["tolerance"] = 1e-8;
//         params["Hypre"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_TRILINOS
// TEST_CASE("Trilinos_b2", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "gr_30_30.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["block_size"] = 2;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_TRILINOS
// TEST_CASE("Trilinos_crystm03", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     std::string MatrixName = "crystm03.mtx";
//     Eigen::SparseMatrix<double> A;
//     loadSymmetric(A, path + "/" + MatrixName);
//     std::cout << "Matrix Load OK" << std::endl;
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["block_size"] = 3;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif

// #ifdef POLYSOLVE_WITH_TRILINOS
// TEST_CASE("Trilinos_smallscale", "[solver]")
// {
//     const std::string path = POLYSOLVE_DATA_DIR;
//     Eigen::SparseMatrix<double> A;
//     const bool ok = loadMarket(A, path + "/A_2.mat");
//     REQUIRE(ok);
//     Eigen::VectorXd b(A.rows());
//     b.setOnes();
//     Eigen::VectorXd x(b.size());
//     x.setZero();
//     Eigen::VectorXd x_b(b.size());
//     x_b.setZero();
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Scalar Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     {
//         clock_t start, end;
//         json solver_info;
//         start = clock();
//         auto solver = LinearSolver::create("Trilinos", "");
//         json params;
//         params["Trilinos"]["block_size"] = 3;
//         params["Trilinos"]["tolerance"] = 1e-8;
//         params["Trilinos"]["max_iter"] = 1000;
//         solver->setParameters(params);
//         solver->analyzePattern(A, A.rows());
//         solver->factorize(A);
//         solver->solve(b, x_b);
//         end = clock();
//         solver->getInfo(solver_info);
//         std::cout << "Block Running time is " << double(end - start) / CLOCKS_PER_SEC << std::endl;
//         std::cout << solver_info["num_iterations"] << std::endl;
//         std::cout << solver_info["final_res_norm"] << std::endl;
//     }
//     const double err = (A * x - b).norm() / b.norm();
//     const double err_b = (A * x_b - b).norm() / b.norm();
//     std::cout << "Scalar relative error " << err << std::endl;
//     std::cout << "Block relative error " << err_b << std::endl;
//     REQUIRE(err < 1e-8);
//     REQUIRE(err_b < 1e-8);
// }
// #endif
