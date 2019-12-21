
#include <catch2/catch.hpp>

#include "../src/stdafx.h"
#include "../src/solnp.hpp"

dlib::matrix<double, 4, 1> powell(const dlib::matrix<double, 5, 1> &m)
/*
This function computes what the Powell function from the
original documentaiton of SOLNP.
*/
{
    /* 4 parameters*/
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);
    const double x4 = m(3);
    const double x5 = m(4);


    // compute the alkyla function and return the result, equality constraint results and the inequality constraint results
    dlib::matrix<double, 4, 1> return_values(4);
    // Function value
    return_values(0) = std::exp(x1 * x2 * x3 * x4 * x5);
    // Equality constraints
    return_values(1) = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 - 10;
    return_values(2) = x2 * x3 - 5 * x4 * x5;
    return_values(3) = x1 * x1 * x1 + x2 * x2 * x2 + 1;
    return return_values;
}


struct powell_functor {
public:
    powell_functor() {};

    dlib::matrix<double, 4, 1> operator()(const dlib::matrix<double, 5, 1> &x) {
        return powell(x);
    }
};

dlib::matrix<double, 8, 1> alkyla(const dlib::matrix<double, 10, 1> &m)
/*
This function computes what the alkyla function from the
original documentaiton of SOLNP.
*/
{
    /* 10 parameters*/
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);
    const double x4 = m(3);
    const double x5 = m(4);
    const double x6 = m(5);
    const double x7 = m(6);
    const double x8 = m(7);
    const double x9 = m(8);
    const double x10 = m(9);


    // compute the alkyla function and return the result, equality constraint results and teh inequality constraint results
    dlib::matrix<double, 8, 1> return_values(8);
    // Function value
    return_values(0) = -0.63 * x4 * x7 + 50.4 * x1 + 3.5 * x2 + x3 + 33.6 * x5;
    // Equality constraints
    return_values(1) = 98.0 * x3 - 0.1 * x4 * x6 * x9 - x3 * x6;
    return_values(2) = 1000.0 * x2 + 100.0 * x5 - 100.0 * x1 * x8;
    return_values(3) = 122.0 * x4 - 100.0 * x1 - 100.0 * x5;
    // Inequality constraints
    return_values(4) = (1.12 * x1 + 0.13167 * x1 * x8 - 0.00667 * x1 * x8 * x8) / x4;
    return_values(5) = (1.098 * x8 - 0.038 * x8 * x8 + 0.325 * x6 + 57.25) / x7;
    return_values(6) = (-0.222 * x10 + 35.82) / x9;
    return_values(7) = (3.0 * x7 - 133.0) / x10;
    return return_values;
}


struct alkyla_functor {
public:
    alkyla_functor() {};

    dlib::matrix<double, 8, 1> operator()(const dlib::matrix<double, 10, 1> &x) {
        return alkyla(x);
    }
};


TEST_CASE("Calculate the Powell function", "[powell]") {
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            -2.0,
            2.0,
            2.0,
            -1,
            -1;

    dlib::matrix<double, 4, 1> result = powell(parameter_data);

    CHECK(result(0) == Approx(0.000335462627903));
    CHECK(result(1) == Approx(4.0));
    CHECK(result(2) == Approx(-1.0));
    CHECK(result(3) == Approx(1.0));


}


TEST_CASE("Optimize the Powell function", "[powell]") {


    /* x0 */
    dlib::matrix<double, 5, 1> parameter_data;
    parameter_data =
            -2.0,
            2.0,
            2.0,
            -1,
            -1;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the equality constraints are all 0
    dlib::matrix<double, 0, 1> constraints = powell(result);
    for (auto row = 1; row < 3; ++row) {
        CHECK(constraints(row) == Approx(0.0));
    }

    // Check the parameters
    CHECK(result(0) == Approx(-1.717142506313303));
    CHECK(result(1) == Approx(1.595708459713134));
    CHECK(result(2) == Approx(1.827247731350245));
    CHECK(result(3) == Approx(-0.763643197991088));
    CHECK(result(4) == Approx(-0.763643197980140));

    REQUIRE(calculate <= 0.053949846871732);

}


TEST_CASE("Calculate the Alkyla function", "[alkyla]") {
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

    CHECK(result(0) == Approx(-59.175999999999931));
    CHECK(result(1) == Approx(4.639999999999418));
    CHECK(result(2) == Approx(14.0));
    CHECK(result(3) == Approx(-58.999999999999773));
    CHECK(result(4) == Approx(1.015869200000000));
    CHECK(result(5) == Approx(0.997758620689655));
    CHECK(result(6) == Approx(0.391666666666666));
    CHECK(result(7) == Approx(0.938064516129032));


}


TEST_CASE("Optimize the Alkyla function", "[alkyla]") {


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

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(alkyla_functor(), parameter_data, ib, logger, 0.0);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    for (auto row = 0; row < parameter_data.nr(); ++row) {
        CHECK(result(row) >= Approx(parameter_data(row, 1)).epsilon(0.01));
        CHECK(result(row) <= Approx(parameter_data(row, 2)).epsilon(0.01));
    }

    dlib::matrix<double, 0, 1> constraints = alkyla(result);
    for (auto row = 1; row < 4; ++row) {
        CHECK(std::abs(constraints(row)) == Approx(0.0).epsilon(0.01));
    }
    for (auto row = 0; row < ib.nr(); ++row) {
        CHECK(constraints(4 + row) >= Approx(ib(row, 0)).epsilon(0.01));
        CHECK(constraints(4 + row) <= Approx(ib(row, 1)).epsilon(0.01));
    }

    // Validate values
    CHECK(result(0) == Approx(16.996376587648808).epsilon(0.01));
    CHECK(result(1) == Approx(15.999402679162111).epsilon(0.01));
    CHECK(result(2) == Approx(57.688358424585260).epsilon(0.01));
    CHECK(result(3) == Approx(30.324890354969355).epsilon(0.01));
    CHECK(result(4) == Approx(19.999989645413820).epsilon(0.01));
    CHECK(result(5) == Approx(90.565424808707550).epsilon(0.01));
    CHECK(result(6) == Approx(94.999992714259020).epsilon(0.01));
    CHECK(result(7) == Approx(10.590140523335723).epsilon(0.01));
    CHECK(result(8) == Approx(1.561646284077410).epsilon(0.01));
    CHECK(result(9) == Approx(1.535353201975077e+02).epsilon(0.01));

    REQUIRE(calculate <= -172.64179);

}