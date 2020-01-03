
#include <catch2/catch.hpp>

#include "../src/solver.cpp"

double alkyla_objective(const std::vector<double> &m)
/*
Alkyla function from the original documentation of SOLNP:
 ALKYLA: this problem has both equality and inequality constraints, as well as the parameter bounds.
*/
{
    const double &x1 = m[0];
    const double &x2 = m[1];
    const double &x3 = m[2];
    const double &x4 = m[3];
    const double &x5 = m[4];
    const double &x6 = m[5];
    const double &x7 = m[6];
    const double &x8 = m[7];
    const double &x9 = m[8];
    const double &x10 = m[9];


    // Function value
    double return_value = -0.63 * x4 * x7 + 50.4 * x1 + 3.5 * x2 + x3 + 33.6 * x5;
    return return_value;
}

std::vector<double> alkyla_equality(const std::vector<double> &m)
/*
Alkyla function from the original documentation of SOLNP:
 ALKYLA: this problem has both equality and inequality constraints, as well as the parameter bounds.
*/
{
    const double &x1 = m[0];
    const double &x2 = m[1];
    const double &x3 = m[2];
    const double &x4 = m[3];
    const double &x5 = m[4];
    const double &x6 = m[5];
    const double &x7 = m[6];
    const double &x8 = m[7];
    const double &x9 = m[8];
    const double &x10 = m[9];


    // compute the alkyla function and return the result, equality constraint results and the inequality constraint results
    std::vector<double> return_values(3);
    // Equality constraints
    return_values[0] = 98.0 * x3 - 0.1 * x4 * x6 * x9 - x3 * x6;
    return_values[1] = 1000.0 * x2 + 100.0 * x5 - 100.0 * x1 * x8;
    return_values[2] = 122.0 * x4 - 100.0 * x1 - 100.0 * x5;
    return return_values;
}


std::vector<double> alkyla_inequality(const std::vector<double> &m)
/*
Alkyla function from the original documentation of SOLNP:
 ALKYLA: this problem has both equality and inequality constraints, as well as the parameter bounds.
*/
{
    const double &x1 = m[0];
    const double &x2 = m[1];
    const double &x3 = m[2];
    const double &x4 = m[3];
    const double &x5 = m[4];
    const double &x6 = m[5];
    const double &x7 = m[6];
    const double &x8 = m[7];
    const double &x9 = m[8];
    const double &x10 = m[9];


    // compute the alkyla function and return the result, equality constraint results and the inequality constraint results
    std::vector<double> return_values(4);
    // Inequality constraints
    return_values[0] = (1.12 * x1 + 0.13167 * x1 * x8 - 0.00667 * x1 * x8 * x8) / x4;
    return_values[1] = (1.098 * x8 - 0.038 * x8 * x8 + 0.325 * x6 + 57.25) / x7;
    return_values[2] = (-0.222 * x10 + 35.82) / x9;
    return_values[3] = (3.0 * x7 - 133.0) / x10;
    return return_values;
}


TEST_CASE("Optimize the Alkyla function pysolnp", "[alkyla]") {


    /* x0, lower, upper */
    std::vector<double> starting_point{17.45, 12.0, 110.0, 30.0, 19.74, 89.2, 92.8, 8.0, 3.6, 155.0};
    std::vector<double> lower_bound{0.0, 0.0, 0.0, 10.0, 0.0, 85.0, 10.0, 3.0, 1.0, 145.0};
    std::vector<double> upper_bound{20.0, 16.0, 120.0, 50.0, 20.0, 93.0, 95.0, 12.0, 4.0, 162.0};

    std::vector<double> eq_value{0, 0, 0};

    /* Inequality function constraints.*/
    std::vector<double> lb{.99, .99, .9, .99};

    std::vector<double> ub{100.0 / 99.0, 100.0 / 99.0, 10.0 / 9.0, 100.0 / 99.0};

    cppsolnp::SolverResult result = solve(
            starting_point,
            lower_bound,
            upper_bound,
            alkyla_objective,
            alkyla_equality,
            eq_value,
            alkyla_inequality,
            lb,
            ub,
            0.0,
            10,
            10,
            1e-5,
            1e-4,
            true);




    // Validate values
    CHECK(result.optimum(0) == Approx(0.169963765876488e2));
    CHECK(result.optimum(1) == Approx(0.159994026791621e2));
    CHECK(result.optimum(2) == Approx(0.576883584245853e2));
    CHECK(result.optimum(3) == Approx(0.303248903549694e2));
    CHECK(result.optimum(4) == Approx(0.199999896454138e2));
    CHECK(result.optimum(5) == Approx(0.905654248087076e2));
    CHECK(result.optimum(6) == Approx(0.949999927142590e2));
    CHECK(result.optimum(7) == Approx(0.105901405233357e2));
    CHECK(result.optimum(8) == Approx(0.015616462840774e2));
    CHECK(result.optimum(9) == Approx(1.535353201975077e2));

    REQUIRE(result.solve_value <= Approx(-1.726412486481025e2));

}
