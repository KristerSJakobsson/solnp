
#include <catch2/catch.hpp>

#include "../src/stdafx.h"
#include "../src/solnp.hpp"

// Benchmark functions
#include "../src/benchmark/alkyla.hpp"
#include "../src/benchmark/box.hpp"
#include "../src/benchmark/entropy.hpp"
#include "../src/benchmark/powell.hpp"
#include "../src/benchmark/wright_four.hpp"
#include "../src/benchmark/wright_nine.hpp"


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
        CHECK(result(row) >= Approx(parameter_data(row, 1)).margin(0.00001));
        CHECK(result(row) <= Approx(parameter_data(row, 2)).margin(0.00001));
    }

    dlib::matrix<double, 0, 1> constraints = alkyla(result);
    for (auto row = 1; row < 4; ++row) {
        CHECK(constraints(row) == Approx(0.0).margin(0.00001));
    }
    for (auto row = 0; row < ib.nr(); ++row) {
        CHECK(constraints(4 + row) >= Approx(ib(row, 0)).margin(0.00001));
        CHECK(constraints(4 + row) <= Approx(ib(row, 1)).margin(0.00001));
    }

    // Validate values
    CHECK(result(0) == Approx(16.996376587648808).margin(0.00001));
    CHECK(result(1) == Approx(15.999402679162111).margin(0.00001));
    CHECK(result(2) == Approx(57.688358424585260).margin(0.00001));
    CHECK(result(3) == Approx(30.324890354969355).margin(0.00001));
    CHECK(result(4) == Approx(19.999989645413820).margin(0.00001));
    CHECK(result(5) == Approx(90.565424808707550).margin(0.00001));
    CHECK(result(6) == Approx(94.999992714259020).margin(0.00001));
    CHECK(result(7) == Approx(10.590140523335723).margin(0.00001));
    CHECK(result(8) == Approx(1.561646284077410).margin(0.00001));
    CHECK(result(9) == Approx(1.535353201975077e+02).margin(0.00001));

    REQUIRE(calculate <= Approx(-172.6412486481025).margin(0.00000001));

}


TEST_CASE("Calculate the Box function", "[box]") {
    dlib::matrix<double, 3, 1> parameter_data;
    parameter_data =
            1.1,
            1.1,
            9.0;

    dlib::matrix<double, 2, 1> result = box(parameter_data);

    CHECK(result(0) == Approx(-10.890000000000002));
    CHECK(result(1) == Approx(-55.560000000000002));


}

TEST_CASE("Optimize the Box function (case a)", "[box]") {


    /* x0 */
    dlib::matrix<double, 3, 3> parameter_data;
    parameter_data =
            1.1, 1.0, 10.0,
            1.1, 1.0, 10.0,
            9.0, 1.0, 10.0;

    dlib::matrix<double, 0, 0> ib;

    cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

    double calculate = cppsolnp::solnp(box_functor(), parameter_data, ib, logger);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the equality constraints are all 0
    dlib::matrix<double, 0, 1> constraints = box(result);
    CHECK(constraints(1) == Approx(0.0).margin(0.000001));

    // Check the parameters
    CHECK(result(0) == Approx(2.8867750701860793));
    CHECK(result(1) == Approx(2.8867750712542191));
    CHECK(result(2) == Approx(5.773407750260735));

    REQUIRE(calculate <= -48.112522068150462);

}


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


TEST_CASE("Optimize the Powell function (rho == 0)", "[powell]") {


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

    double calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 0.0);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the equality constraints are all 0
    dlib::matrix<double, 0, 1> constraints = powell(result);
    for (auto row = 1; row < 3; ++row) {
        CHECK(constraints(row) == Approx(0.0).margin(0.000001));
    }

    // Check the parameters
    CHECK(result(0) == Approx(-1.717142506313303));
    CHECK(result(1) == Approx(1.595708459713134));
    CHECK(result(2) == Approx(1.827247731350245).margin(0.0001));
    CHECK(result(3) == Approx(-0.763643197991088));
    CHECK(result(4) == Approx(-0.763643197980140));

    REQUIRE(calculate <= 0.053949846871732);

}

TEST_CASE("Optimize the Powell function (rho == 1)", "[powell]") {


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

    double calculate = cppsolnp::solnp(powell_functor(), parameter_data, ib, logger, 1.0);

    dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

    // Check the equality constraints are all 0
    dlib::matrix<double, 0, 1> constraints = powell(result);
    for (auto row = 1; row < 3; ++row) {
        CHECK(constraints(row) == Approx(0.0).margin(0.000001));
    }

    // Check the parameters
    CHECK(result(0) == Approx(-1.717142506313303));
    CHECK(result(1) == Approx(1.595708459713134));
    CHECK(result(2) == Approx(1.827247731350245));
    CHECK(result(3) == Approx(-0.763643197991088));
    CHECK(result(4) == Approx(-0.763643197980140));

    REQUIRE(calculate <= 0.053949846871732);

}
