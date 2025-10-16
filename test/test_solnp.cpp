#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "stdafx.h"
#include "solnp.hpp"

// Benchmark functions
#include "./benchmark/alkyla.hpp"
#include "./benchmark/box.hpp"
#include "./benchmark/entropy.hpp"
#include "./benchmark/powell.hpp"
#include "./benchmark/rosen_suzuki.hpp"
#include "./benchmark/wright_four.hpp"
#include "./benchmark/wright_nine.hpp"
#include "catch2/matchers/catch_matchers.hpp"


TEST_CASE("Calculate trivial Quadratic function", "[y=x^2][solnp][sanity]")
{
    dlib::matrix<double, 1, 1> parameter_data;
    parameter_data = 1.0;

    auto quadratic_function = [](const dlib::matrix<double, 1, 1>& m) -> dlib::matrix<double, 1, 1>
    {
        return dlib::mat(m(0) * m(0));
    };

    cppsolnp::SolveResult calculate = cppsolnp::solnp(quadratic_function, parameter_data);

    CHECK(calculate.converged == true);

    // Validate values
    dlib::matrix<double, 0, 1> result = calculate.optimum;
    CHECK(result(0) == Catch::Approx(0.0).margin(1e-2));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.0).margin(1e-3));
}

TEST_CASE("Throws when tolerance is not finite", "[solnp][exception]") {
    dlib::matrix<double, 2, 1> param;
    param = 1.0, 2.0;
    auto functor = [](const dlib::matrix<double>& x) {
        dlib::matrix<double, 2, 1> y;
        return y;
    };
    double inf_tol = std::numeric_limits<double>::infinity();

    REQUIRE_THROWS_AS(
        cppsolnp::solnp(functor, param, nullptr, 1.0, 400, 800, 1e-7, inf_tol),
        std::invalid_argument
    );
}

TEST_CASE("Throws exception for contradicting equality constraints", "[y=x^2][solnp][exception]")
{
    dlib::matrix<double, 1, 1> parameter_data;
    parameter_data = 1.0;

    auto quadratic_function = [](const dlib::matrix<double, 1, 1>& m) -> dlib::matrix<double, 3, 1>
    {
        const double x1 = m(0);

        // compute the box function and return the result, equality constraint results and the equality constraint results
        dlib::matrix<double, 3, 1> return_values(3);
        // Function value
        return_values(0) = x1 * x1; // x1^2
        // Equality constraints
        return_values(1) = x1 - 1.0; // x1=1
        return_values(2) = x1 - 2.0; // x1=2
        return return_values;
    };

    REQUIRE_THROWS_WITH(cppsolnp::solnp(quadratic_function, parameter_data),
                        "Encountered Singular matrix when trying to solve. This can happen for example if you have contradicting equality constraints.")
    ;
}

// TEST_CASE("Throws exception for zero tolerance", "[solnp][exception]")
// {
//     dlib::matrix<double, 1, 1> parameter_data;
//     parameter_data = 1.0;
//     dlib::matrix<double> hessian_matrix = dlib::zeros_matrix<double>(1, 4);
//     auto functor = [](const dlib::matrix<double, 1, 1>& m) -> dlib::matrix<double, 2, 1>
//     {
//         dlib::matrix<double, 2, 1> return_values(2);
//         return_values(0) = dlib::mat(m(0) * m(0));
//         return_values(1) = m(0) * 5 - 5;
//         return return_values;
//     };
//     std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();
//     REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data,  logger, 1.0, 400, 800, 1E-07, 0),
//                         "Encountered Singular matrix when trying to solve. This can happen for example if you have contradicting equality constraints.");
//
// }

TEST_CASE("Throws exception for too many parameter columns", "[solnp][exception]")
{
    dlib::matrix<double> parameter_data(2, 4); // Intentionally set size on runtime to bypass template checks
    parameter_data = 1, 2, 3, 4, 5, 6, 7, 8;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "Parameter array must have three columns or less.");
}

TEST_CASE("Throws exception for invalid parameter bounds (2 columns)", "[solnp][exception]")
{
    dlib::matrix<double, 2, 2> parameter_data;
    parameter_data = 1.0, 1.0, 2.0, 0.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "The lower bounds of the parameter constraints must be strictly less than the upper bounds.");
}

TEST_CASE("Throws exception for invalid parameter bounds (3 columns)", "[solnp][exception]")
{
    dlib::matrix<double, 2, 3> parameter_data;
    parameter_data = 1.0, 2.0, 2.0, 2.0, 2.0, 1.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "The lower bounds of the parameter constraints must be strictly less than the upper bounds.");
}

TEST_CASE("Throws exception for initial parameter out of bounds (3 columns)", "[solnp][exception]")
{
    dlib::matrix<double, 2, 3> parameter_data;
    parameter_data = 5.0, 0.0, 2.0, 5.0, 0.0, 2.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "Initial parameter values must be within the bounds.");

    parameter_data = -1.0, 0.0, 2.0, -1.0, 0.0, 2.0;
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "Initial parameter values must be within the bounds.");
}

TEST_CASE("Throws exception for initial inequalities out of bounds", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    dlib::matrix<double> ib(2, 3); // Intentionally set size on runtime to bypass template checks
    ib = 5.0, 0.0, 2.0, 5.0, 0.0, 2.0;
    auto functor = [](const dlib::matrix<double>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib),
                        "Initial inequalities must be within bounds.");

    ib = -1.0, 0.0, 2.0, -1.0, 0.0, 2.0;
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib),
                        "Initial inequalities must be within bounds.");
}

TEST_CASE("Throws exception for invalid inequality bounds (2 columns)", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    dlib::matrix<double, 2, 2> ib;
    ib = 1.0, 1.0, 2.0, 0.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib),
                        "The lower bounds of the inequality constraints must be strictly less than the upper bounds.");
}

TEST_CASE("Throws exception for too many inequality columns", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    dlib::matrix<double> ib(2, 4); // Intentionally set size on runtime to bypass template checks
    ib = 1, 2, 3, 4, 5, 6, 7, 8;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib),
                        "Inequality constraints must have 2 or 3 columns.");
}

TEST_CASE("Throws exception for invalid hessian matrix size", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    dlib::matrix<double, 2, 2> ib;
    ib = 0.0, 1.0, 0.0, 1.0;
    dlib::matrix<double> hessian_matrix = dlib::zeros_matrix<double>(1, 4);
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::mat(m(0) + m(1)); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib, hessian_matrix),
                        "The provided hessian matrix override was of invalid dimension.");

    dlib::matrix<double> hessian_matrix2 = dlib::zeros_matrix<double>(4, 1);
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib, hessian_matrix2),
                        "The provided hessian matrix override was of invalid dimension.");
}

TEST_CASE("Throws exception for cost function with more than one column", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::ones_matrix<double>(1, 2); };
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data),
                        "sqp_min cost function must return only 1 column.");
}

TEST_CASE("Throws exception for cost function with too few constraints", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return dlib::ones_matrix<double>(1, 1); };
    dlib::matrix<double, 2, 2> ib;
    ib = 0.0, 1.0, 0.0, 1.0;
    REQUIRE_THROWS_WITH(cppsolnp::solnp(functor, parameter_data, ib),
                        "sqp_min the number of constraints in the cost function does not match the call to sqp_min.");
}

TEST_CASE("Fails to reach convergence", "[solnp][exception]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 1.0, 2.0;
    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();
    auto functor = [](const dlib::matrix<double, 2, 1>& m) { return m; };
    cppsolnp::SolveResult calculate =  cppsolnp::solnp(functor, parameter_data, logger, 1.0, 0, 0, 1e-7, 0);

    // Expect this to not converge since 0 iterations are allowed and tolerance is 0
    CHECK(calculate.converged == false);
}

TEST_CASE("Fails gracefully for contradicting inequality constraints", "[y=x^2]")
{
    dlib::matrix<double, 1, 1> parameter_data;
    parameter_data = 1.0;

    auto quadratic_function = [](const dlib::matrix<double, 1, 1>& m) -> dlib::matrix<double, 3, 1>
    {
        const double x1 = m(0);

        // compute the box function and return the result, equality constraint results and the equality constraint results
        dlib::matrix<double, 3, 1> return_values(3);
        // Function value
        return_values(0) = x1 * x1; // x1^2
        // Inequality constraints
        return_values(1) = x1; // x1=1
        return_values(2) = x1; // x1=2
        return return_values;
    };

    dlib::matrix<double, 2, 2> ib;
    ib =
        0, 1,
        3, 4;

    cppsolnp::SolveResult calculate = cppsolnp::solnp(quadratic_function, parameter_data, ib);

    CHECK(calculate.converged == false);
}

TEST_CASE("Solves function with lower and upper bounds", "[f(x,y) = (x-1.5)^2  + (y-1.5)^2]")
{
    dlib::matrix<double, 2, 2> parameter_data;
    parameter_data =
        0.0, 1.0,
        1.0, 2.0;

    auto functor = [](const dlib::matrix<double, 2, 1>& m) -> dlib::matrix<double, 1, 1>
    {
        dlib::matrix<double, 1, 1> out;
        out(0) = std::pow(m(0) - 1.5, 2) + std::pow(m(1) - 1.5, 2);
        return out;
    };

    auto result = cppsolnp::solnp(functor, parameter_data);

    /* Function f(x,y) = (x-1.5)^2  + (y-1.5)^2 has minimum in x=1.5, y=1.5,
       however bounds limit x to 1.0, giving minimum in x=1.0, y=1.5 */
    CHECK(result.converged == true);
    CHECK(result.optimum(0) == Catch::Approx(1.0));
    CHECK(result.optimum(1) == Catch::Approx(1.5));
}

TEST_CASE("One-column inequality matrix is treated as no inequalities", "[solnp][inequality][edge]")
{
    dlib::matrix<double, 2, 1> parameter_data;
    parameter_data = 3.0, -1.0;

    // One-column inequality matrix (should be treated as no inequalities)
    dlib::matrix<double> ib(2, 1);
    ib = 0.5, 1.5;

    // Provide explicit hessian matching expected size when inequalities are ignored (2 parameters + 0 inequalities)
    dlib::matrix<double> hessian = dlib::identity_matrix<double>(parameter_data.nr());

    auto functor = [](const dlib::matrix<double, 2, 1>& m) -> dlib::matrix<double, 1, 1>
    {
        dlib::matrix<double, 1, 1> out;
        out(0) = std::pow(m(0) - 1.0, 2) + std::pow(m(1) - 2.0, 2);
        return out;
    };

    // Result without providing inequality matrix
    cppsolnp::SolveResult res_no_ib = cppsolnp::solnp(functor, parameter_data);

    // Result with one-column inequality matrix and explicit hessian
    cppsolnp::SolveResult res_with_ib = cppsolnp::solnp(functor, parameter_data, ib, hessian);

    // Both should behave the same (inequality column treated as no inequalities)
    CHECK(res_with_ib.converged == res_no_ib.converged);
    REQUIRE(res_with_ib.optimum.nr() == res_no_ib.optimum.nr());
    for (long i = 0; i < res_no_ib.optimum.nr(); ++i) {
        CHECK(res_with_ib.optimum(i) == Catch::Approx(res_no_ib.optimum(i)).margin(1e-6));
    }
    CHECK(res_with_ib.solve_value == Catch::Approx(res_no_ib.solve_value).margin(1e-6));
}

/*
 * Below tests compare the scenarios and results from the original SOLNP in Matlab with the results from C++ SOLNP
 * Notably, these differ a bit from what pysolnp will generate due to:
 * - Objective/Constraint functions in pysolnp are subject to rounding errors in Python
 * - The precompiled binaries for pysolnp does not utilize BLAS or LAPACK.
 * Note that depending on if your system has BLAS and/or LAPACK installed the results from these tests might very.
 * */

TEST_CASE("Calculate the Alkyla function", "[alkyla]")
{
    dlib::matrix<double, 10, 1> parameter_data;
    parameter_data =
        17.45,
        12.0,
        110.0,
        30.0,
        19.74,
        89.2,
        92.8,
        8.0,
        3.6,
        155.0;

    dlib::matrix<double, 8, 1> result = alkyla(parameter_data);

    CHECK(result(0) == Catch::Approx(-59.175999999999931));
    CHECK(result(1) == Catch::Approx(4.639999999999418));
    CHECK(result(2) == Catch::Approx(14.0));
    CHECK(result(3) == Catch::Approx(-58.999999999999773));
    CHECK(result(4) == Catch::Approx(1.015869200000000));
    CHECK(result(5) == Catch::Approx(0.997758620689655));
    CHECK(result(6) == Catch::Approx(0.391666666666666));
    CHECK(result(7) == Catch::Approx(0.938064516129032));
}


TEST_CASE("Optimize the Alkyla function manual hessian", "[alkyla]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 10, 3> parameter_data;
    parameter_data =
        17.45, 0.0, 20.0,
        12.0, 0.0, 16.0,
        110.0, 0.0, 120.0,
        30.0, 10.0, 50.0,
        19.74, 0.0, 20.0,
        89.2, 85.0, 93.0,
        92.8, 10.0, 95.0,
        8.0, 3.0, 12.0,
        3.6, 1.0, 4.0,
        155.0, 145.0, 162.0;

    /* Inequality function constraints.*/
    dlib::matrix<double, 4, 2> ib;
    ib =
        .99, 100.0 / 99.0,
        .99, 100.0 / 99.0,
        .9, 10.0 / 9.0,
        .99, 100.0 / 99.0;

    dlib::matrix<double> hessian_matrix;
    hessian_matrix = dlib::identity_matrix<double>(parameter_data.nr() + ib.nr());

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(alkyla_functor(), parameter_data, ib, hessian_matrix, logger, 0.0,
                                                      10, 10, 1e-5,
                                                      1e-4);

    CHECK(calculate.converged == true);

    // Validate values
    dlib::matrix<double, 0, 1> result = calculate.optimum;
    CHECK(result(0) == Catch::Approx(0.169963765876488e2).margin(1e-2));
    CHECK(result(1) == Catch::Approx(0.159994026791621e2).margin(1e-2));
    CHECK(result(2) == Catch::Approx(0.576883584245853e2).margin(1e-2));
    CHECK(result(3) == Catch::Approx(0.303248903549694e2).margin(1e-2));
    CHECK(result(4) == Catch::Approx(0.199999896454138e2).margin(1e-2));
    CHECK(result(5) == Catch::Approx(0.905654248087076e2).margin(1e-2));
    CHECK(result(6) == Catch::Approx(0.949999927142590e2).margin(1e-2));
    CHECK(result(7) == Catch::Approx(0.105901405233357e2).margin(1e-2));
    CHECK(result(8) == Catch::Approx(0.015616462840774e2).margin(1e-2));
    CHECK(result(9) == Catch::Approx(1.535353201975077e2).margin(1e-2));

    REQUIRE(calculate.solve_value <= Catch::Approx(-1.726412486481025e2).margin(1e-3));
}


TEST_CASE("Calculate the Box function", "[box]")
{
    dlib::matrix<double, 3, 1> parameter_data;
    parameter_data =
        1.1,
        1.1,
        9.0;

    dlib::matrix<double, 2, 1> result = box(parameter_data);

    CHECK(result(0) == Catch::Approx(-10.890000000000002));
    CHECK(result(1) == Catch::Approx(-55.560000000000002));
}

TEST_CASE("Optimize the Box function (case a)", "[box]")
{
    /* x0 */
    dlib::matrix<double, 3, 3> parameter_data;
    parameter_data =
        1.1, 1.0, 10.0,
        1.1, 1.0, 10.0,
        9.0, 1.0, 10.0;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(box_functor(), parameter_data, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Check the parameters
    CHECK(result(0) == Catch::Approx(2.886775069536727));
    CHECK(result(1) == Catch::Approx(2.886775072009683));
    CHECK(result(2) == Catch::Approx(5.773407750048355).margin(1e-2));

    REQUIRE(calculate.solve_value <= Catch::Approx(-48.1125220681));
}


TEST_CASE("Optimize the Box function (case b)", "[box]")
{
    /* x0 */
    dlib::matrix<double, 3, 3> parameter_data;
    parameter_data =
        5.5, 1.0, 10.0,
        5.5, 1.0, 10.0,
        5.5, 1.0, 10.0;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(box_functor(), parameter_data, logger, 1.0, 10, 10, 1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Check the parameters
    CHECK(result(0) == Catch::Approx(2.888765743268910).margin(1e-2));
    CHECK(result(1) == Catch::Approx(2.888765747765645).margin(1e-2));
    CHECK(result(2) == Catch::Approx(5.765448483893261).margin(1e-2));

    REQUIRE(calculate.solve_value <= Catch::Approx(-48.112480408240664));
}


TEST_CASE("Calculate the Entropy function", "[entropy]")
{
    dlib::matrix<double, 10, 1> parameter_data;
    parameter_data = 0.8474, 0.4524, 0.8075, 0.4832, 0.6135, 0.2749, 0.8807, 0.6538, 0.4899, 0.7741;

    dlib::matrix<double, 2, 1> result = entropy(parameter_data);

    CHECK(result(0) == Catch::Approx(4.849345605));
    CHECK(result(1) == Catch::Approx(-3.7226));
}


TEST_CASE("Optimize the Entropy function", "[entropy]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 10, 3> parameter_data;
    parameter_data =
        0.8474, 0.0, 10.0,
        0.4524, 0.0, 10.0,
        0.8075, 0.0, 10.0,
        0.4832, 0.0, 10.0,
        0.6135, 0.0, 10.0,
        0.2749, 0.0, 10.0,
        0.8807, 0.0, 10.0,
        0.6538, 0.0, 10.0,
        0.4899, 0.0, 10.0,
        0.7741, 0.0, 10.0;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(entropy_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5,
                                                      1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(0.857717373389226));
    CHECK(result(1) == Catch::Approx(0.857627630428985));
    CHECK(result(2) == Catch::Approx(0.858758912629577));
    CHECK(result(3) == Catch::Approx(0.857611377467009));
    CHECK(result(4) == Catch::Approx(0.857445851804114));
    CHECK(result(5) == Catch::Approx(0.857975885723616));
    CHECK(result(6) == Catch::Approx(2.279903520424556));
    CHECK(result(7) == Catch::Approx(0.857328700404536));
    CHECK(result(8) == Catch::Approx(0.857608122266489));
    CHECK(result(9) == Catch::Approx(0.858022625546437));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.185478885901993));
}


TEST_CASE("Calculate the Powell function", "[powell]")
{
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        -2.0,
        2.0,
        2.0,
        -1,
        -1;

    dlib::matrix<double, 4, 1> result = powell(parameter_data);

    CHECK(result(0) == Catch::Approx(0.000335462627903));
    CHECK(result(1) == Catch::Approx(4.0));
    CHECK(result(2) == Catch::Approx(-1.0));
    CHECK(result(3) == Catch::Approx(1.0));
}

TEST_CASE("Optimize the Powell function (rho == 0)", "[powell]")
{
    /* x0 */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        -2.0,
        2.0,
        2.0,
        -1,
        -1;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 0.0, 10, 10, 1e-5,
                                                      1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Check the parameters
    CHECK(result(0) == Catch::Approx(-1.717126203723513));
    CHECK(result(1) == Catch::Approx(1.595689596511580));
    CHECK(result(2) == Catch::Approx(1.827278075550860));
    CHECK(result(3) == Catch::Approx(-0.763645042210886));
    CHECK(result(4) == Catch::Approx(-0.763645042234952));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.053949827793391));
}

TEST_CASE("Optimize the Powell function (rho == 1)", "[powell]")
{
    /* x0 */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        -2.0,
        2.0,
        2.0,
        -1,
        -1;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5,
                                                      1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Check the parameters
    CHECK(result(0) == Catch::Approx(-1.717142506313303));
    CHECK(result(1) == Catch::Approx(1.595708459713134));
    CHECK(result(2) == Catch::Approx(1.827247731350245));
    CHECK(result(3) == Catch::Approx(-0.763643197991088));
    CHECK(result(4) == Catch::Approx(-0.763643197980140));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.0539498469));
}


TEST_CASE("Calculate the Wright4 function", "[wright4]")
{
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1.0, 3.0, -0.5, -2.0, -3.0;

    dlib::matrix<double, 4, 1> result = wright_four(parameter_data);


    CHECK(result(0) == Catch::Approx(68.9375));
    CHECK(result(1) == Catch::Approx(1.632359313));
    CHECK(result(2) == Catch::Approx(-0.07842712475));
    CHECK(result(3) == Catch::Approx(1.0));
}


TEST_CASE("Optimize the Wright4 function (case a, rho==10)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116639595144975));
    CHECK(result(1) == Catch::Approx(1.220448420527845));
    CHECK(result(2) == Catch::Approx(1.537782093973876));
    CHECK(result(3) == Catch::Approx(1.972752470314671));
    CHECK(result(4) == Catch::Approx(1.791088179957703));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310831271171));
}


TEST_CASE("Optimize the Wright4 function (case a, rho==1)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116643803185402));
    CHECK(result(1) == Catch::Approx(1.220453096801827));
    CHECK(result(2) == Catch::Approx(1.537779890421572));
    CHECK(result(3) == Catch::Approx(1.972741018934506));
    CHECK(result(4) == Catch::Approx(1.791081456246007));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310831002758));
}


TEST_CASE("Optimize the Wright4 function (case a, rho==0)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 1, 1, 1, 1, 1;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116609971954848));
    CHECK(result(1) == Catch::Approx(1.220440489315363));
    CHECK(result(2) == Catch::Approx(1.537788998574534));
    CHECK(result(3) == Catch::Approx(1.972781636643539));
    CHECK(result(4) == Catch::Approx(1.791135716667731));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310831942731));
}

TEST_CASE("Optimize the Wright4 function (case b, rho==10)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116634078024861));
    CHECK(result(1) == Catch::Approx(1.220444234062229));
    CHECK(result(2) == Catch::Approx(1.537784320405761));
    CHECK(result(3) == Catch::Approx(1.972763493598540));
    CHECK(result(4) == Catch::Approx(1.791097035073237));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310831022048));
}


TEST_CASE("Optimize the Wright4 function (case b, rho==1)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116634870722612));
    CHECK(result(1) == Catch::Approx(1.220440576122608));
    CHECK(result(2) == Catch::Approx(1.537785458091246));
    CHECK(result(3) == Catch::Approx(1.972770662225987));
    CHECK(result(4) == Catch::Approx(1.791095775742904));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310830648204));
}


TEST_CASE("Optimize the Wright4 function (case b, rho==0)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = 2, 2, 2, 2, 2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.116636615308806));
    CHECK(result(1) == Catch::Approx(1.220440910051733));
    CHECK(result(2) == Catch::Approx(1.537785097114965));
    CHECK(result(3) == Catch::Approx(1.972769218398222));
    CHECK(result(4) == Catch::Approx(1.791092977407627));

    REQUIRE(calculate.solve_value <= Catch::Approx(0.029310830686406));
}

TEST_CASE("Optimize the Wright4 function (case c, rho==10)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-0.703068933803915));
    CHECK(result(1) == Catch::Approx(2.635653020521741));
    CHECK(result(2) == Catch::Approx(-0.099089129462807));
    CHECK(result(3) == Catch::Approx(-1.797457648464959));
    CHECK(result(4) == Catch::Approx(-2.844671274183172));

    REQUIRE(calculate.solve_value <= Catch::Approx(44.022877145171257));
}


TEST_CASE("Optimize the Wright4 function (case c, rho==1)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-0.703523524834065));
    CHECK(result(1) == Catch::Approx(2.635727880354312));
    CHECK(result(2) == Catch::Approx(-0.096451726879527));
    CHECK(result(3) == Catch::Approx(-1.797997954636873));
    CHECK(result(4) == Catch::Approx(-2.842832966138447));

    REQUIRE(calculate.solve_value <= Catch::Approx(44.022075061138295));
}


TEST_CASE("Optimize the Wright4 function (case c, rho==0)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 3, -0.5, -2, -3;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-0.703376376169597));
    CHECK(result(1) == Catch::Approx(2.635700821034921));
    CHECK(result(2) == Catch::Approx(-0.096657075140750));
    CHECK(result(3) == Catch::Approx(-1.797935433028946));
    CHECK(result(4) == Catch::Approx(-2.843427812167240));

    REQUIRE(calculate.solve_value <= Catch::Approx(44.022128023467303));
}


TEST_CASE("Optimize the Wright4 function (case d, rho==10)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 10, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-1.273052029422237));
    CHECK(result(1) == Catch::Approx(2.410353869751445));
    CHECK(result(2) == Catch::Approx(1.194859244564641));
    CHECK(result(3) == Catch::Approx(-0.154238130707192));
    CHECK(result(4) == Catch::Approx(-1.571027698602370));

    REQUIRE(calculate.solve_value <= Catch::Approx(27.871905223431018));
}


TEST_CASE("Optimize the Wright4 function (case d, rho==1)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-1.272710813834789));
    CHECK(result(1) == Catch::Approx(2.410353128557497));
    CHECK(result(2) == Catch::Approx(1.194780438668736));
    CHECK(result(3) == Catch::Approx(-0.154425728750688));
    CHECK(result(4) == Catch::Approx(-1.571448524460437));

    REQUIRE(calculate.solve_value <= Catch::Approx(27.871903584038883));
}


TEST_CASE("Optimize the Wright4 function (case d, rho==0)", "[wright4]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data = -1, 2, 1, -2, -2;

    dlib::matrix<double, 0, 0> ib;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_four_functor(), parameter_data, ib, logger, 0.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-1.273023411221166));
    CHECK(result(1) == Catch::Approx(2.410362154910699));
    CHECK(result(2) == Catch::Approx(1.194843255003828));
    CHECK(result(3) == Catch::Approx(-0.154284646055543));
    CHECK(result(4) == Catch::Approx(-1.571062971404984));

    REQUIRE(calculate.solve_value <= Catch::Approx(27.871904866028800));
}


TEST_CASE("Calculate the Wright9 function", "[wright9]")
{
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        1.091,
        -3.174,
        1.214,
        -1.614,
        2.134;

    dlib::matrix<double, 4, 1> result = wright_nine(parameter_data);

    CHECK(result(0) == Catch::Approx(-1.815401107e+3));
    CHECK(result(1) == Catch::Approx(0.019897305e+3));
    CHECK(result(2) == Catch::Approx(-0.001999274866e+3));
    CHECK(result(3) == Catch::Approx(0.007022058536e+3));
}


TEST_CASE("Optimize the Wright9 function (case a, rho==1)", "[wright9]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        1,
        1,
        1,
        1,
        1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
        -100, 20,
        -2, 100,
        5, 100;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-0.082010879422399));
    CHECK(result(1) == Catch::Approx(3.692422439415791));
    CHECK(result(2) == Catch::Approx(2.487343192052682));
    CHECK(result(3) == Catch::Approx(0.377176847262690));
    CHECK(result(4) == Catch::Approx(0.173650632156765));

    REQUIRE(calculate.solve_value <= Catch::Approx(-2.104078394423900e2));
}


TEST_CASE("Optimize the Wright9 function (case a, rho==100)", "[wright9]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        1,
        1,
        1,
        1,
        1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
        -100, 20,
        -2, 100,
        5, 100;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 100.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(-0.081246392868120));
    CHECK(result(1) == Catch::Approx(3.689847468079701));
    CHECK(result(2) == Catch::Approx(2.491144012826899));
    CHECK(result(3) == Catch::Approx(0.377589784583978));
    CHECK(result(4) == Catch::Approx(0.173399403653944));

    REQUIRE(calculate.solve_value <= Catch::Approx(-2.104073432066561e2));
}


TEST_CASE("Optimize the Wright9 function (case b, rho==100)", "[wright9]")
{
    /* x0, lower, upper */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
        1.091,
        -3.174,
        1.214,
        -1.614,
        2.134;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
        -100, 20,
        -2, 100,
        5, 100;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(wright_nine_functor(), parameter_data, ib, logger, 100.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(1.479634676414054));
    CHECK(result(1) == Catch::Approx(-2.636607129671721));
    CHECK(result(2) == Catch::Approx(1.054666762278611));
    CHECK(result(3) == Catch::Approx(-1.611508943269783));
    CHECK(result(4) == Catch::Approx(2.673892424752704));

    REQUIRE(calculate.solve_value <= Catch::Approx(-2.500584227790517e3));
}


TEST_CASE("Calculate the Rosen-Suzuki function", "[rosen_suzuki]")
{
    dlib::matrix<double, 4, 1> parameter_data;
    parameter_data = 1, 1, 1, 1;

    dlib::matrix<double, 4, 1> result = rosen_suzuki(parameter_data);

    CHECK(result(0) == Catch::Approx(-19.0));
    CHECK(result(1) == Catch::Approx(4.0));
    CHECK(result(2) == Catch::Approx(6.0));
    CHECK(result(3) == Catch::Approx(1.0));
}

TEST_CASE("Optimize the Rosen-Suzuki function", "[rosen_suzuki]")
{
    /* x0 */
    dlib::matrix<double, 4, 1> parameter_data;
    parameter_data = 1, 1, 1, 1;

    /* Inequality function constraints.*/
    dlib::matrix<double, 3, 2> ib;
    ib =
        0, 1000,
        0, 1000,
        0, 1000;

    std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

    cppsolnp::SolveResult calculate = cppsolnp::solnp(rosen_suzuki_functor(), parameter_data, ib, logger, 1.0, 10, 10,
                                                      1e-5, 1e-4);

    dlib::matrix<double, 0, 1> result = calculate.optimum;

    // Validate values
    CHECK(result(0) == Catch::Approx(0.000230374253029));
    CHECK(result(1) == Catch::Approx(0.998564179709879));
    CHECK(result(2) == Catch::Approx(2.000278419459943));
    CHECK(result(3) == Catch::Approx(-0.999859532645383));

    REQUIRE(calculate.solve_value <= Catch::Approx(-43.999759237182886));
}
