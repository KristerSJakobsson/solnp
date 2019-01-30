#include "../src/stdafx.h"
#include <catch2/catch.hpp>

#include "../src/solnp.hpp"

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


    // compute the alkyla function and return the result, equality restraint results and teh inequality restraint results
    dlib::matrix<double, 8, 1> return_values(8);
    // Function value
    return_values(0) = -0.63 * x4 * x7 + 50.4 * x1 + 3.5 * x2 + x3 + 33.6 * x5;
    // Equality constraints
    return_values(1) = 98 * x3 - 0.1 * x4 * x6 * x9 - x3 * x6;
    return_values(2) = 1000 * x2 + 100 * x5 - 100 * x1 * x8;
    return_values(3) = 122 * x4 - 100 * x1 - 100 * x5;
    // Inequality constraints
    return_values(4) = (1.12 * x1 + 0.13167 * x1 * x8 - 0.00667 * x1 * x8 * x8) / x4;
    return_values(5) = (1.098 * x8 - 0.038 * x8 * x8 + 0.325 * x6 + 57.25) / x7;
    return_values(6) = (-0.222 * x10 + 35.82) / x9;
    return_values(7) = (3 * x7 - 133) / x10;
    return return_values;
}


struct alkyla_functor {
public:
    alkyla_functor() {};

    dlib::matrix<double, 8, 1> operator()(const dlib::matrix<double, 10, 1> &x) {
        return alkyla(x);
    }
};


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
            155, 145.0, 162.0;

    /* Inequality function constraints.*/
    dlib::matrix<double, 4, 2> ib;
    ib =
            .99, 100.0 / 99.0,
            .99, 100.0 / 99.0,
            .9, 10.0 / 9.0,
            .99, 100.0 / 99.0;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(alkyla_functor(), parameter_data, ib, logger);

    REQUIRE(calculate <= -172.64179);

}