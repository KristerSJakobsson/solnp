
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "stdafx.h"
#include "utils.hpp"


TEST_CASE("Conditional number statically sized", "[utils]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 1.0, 3.0,
            2.0, 6.0, 8.0,
            6.0, 8.0, 18.0;

    double cond =
            cppsolnp::conditional_number(A);

    CHECK(cond == Catch::Approx(31.413545501499051));

}

TEST_CASE("Conditional number dynamically sized", "[utils]") {
    dlib::matrix<double> A(3, 3);
    A =
            21.0, 14.0, -2.0,
            2.0, 1.0, -3.0,
            -21.0, 0.0, 18.0;

    double cond =
            cppsolnp::conditional_number(A);

    CHECK(cond == Catch::Approx(25.003922978899976));

}

TEST_CASE("Euclidean norm statically sized", "[utils]") {
    dlib::matrix<double, 4, 1> x;
    x =
            -2.0, 2.0, 3.0, 3.5;

    double norm = cppsolnp::euclidean_norm(x);

    CHECK(norm == Catch::Approx(5.408326913195984));

}

TEST_CASE("Euclidean norm dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(4);
    x =
            9.0, 15.0, 3.0, -2.5;

    double norm = cppsolnp::euclidean_norm(x);

    CHECK(norm == Catch::Approx(17.923448328934921));

}

TEST_CASE("Infinity norm statically sized", "[utils]") {
    dlib::matrix<double, 4, 1> x;
    x =
            -5.5, 2.0, 3.0, 3.5;

    double norm = cppsolnp::infinity_norm(x);

    CHECK(norm == Catch::Approx(5.5));

}

TEST_CASE("Infinity norm dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(4);
    x =
            9.0, 15.0, 3.0, -2.5;

    double norm = cppsolnp::infinity_norm(x);

    CHECK(norm == Catch::Approx(15.0));

}

TEST_CASE("Pointwise divide statically sized", "[utils]") {
    dlib::matrix<double, 3, 1> x;
    x =
            2.0, 1.0, 3.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 3, 1> x_divide_y;
    x_divide_y =
            cppsolnp::pointwise_divide(x, y);

    CHECK(x_divide_y(0) == Catch::Approx(x(0) / y(0)));
    CHECK(x_divide_y(1) == Catch::Approx(x(1) / y(1)));
    CHECK(x_divide_y(2) == Catch::Approx(x(2) / y(2)));

}

TEST_CASE("Pointwise divide dynamically sized", "[utils]") {
    dlib::matrix<double, 3, 1> x;
    x =
            2.0, 1.0, 3.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 3, 1> x_divide_y;
    x_divide_y =
            cppsolnp::pointwise_divide(x, y);

    CHECK(x_divide_y(0) == Catch::Approx(x(0) / y(0)));
    CHECK(x_divide_y(1) == Catch::Approx(x(1) / y(1)));
    CHECK(x_divide_y(2) == Catch::Approx(x(2) / y(2)));

}

TEST_CASE("Elementwise max statically sized", "[utils]") {
    dlib::matrix<double, 6, 1> x;
    x =
            2.0, -1.0, 3.0, 5.0, 0.0, 1.5;

    dlib::matrix<double, 6, 1> y;
    y =
            cppsolnp::elementwise_max(x, 2.0);

    CHECK(y(0) == Catch::Approx(std::max(x(0), 2.0)));
    CHECK(y(1) == Catch::Approx(std::max(x(1), 2.0)));
    CHECK(y(2) == Catch::Approx(std::max(x(2), 2.0)));
    CHECK(y(3) == Catch::Approx(std::max(x(3), 2.0)));
    CHECK(y(4) == Catch::Approx(std::max(x(4), 2.0)));
    CHECK(y(5) == Catch::Approx(std::max(x(5), 2.0)));

}

TEST_CASE("Elementwise max dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(6);
    x =
            1.5, 2.0, -99.0, 15.0, 0.15, 1.5;

    dlib::matrix<double, 0, 1> y;
    y =
            cppsolnp::elementwise_max(x, 2.0);

    CHECK(y(0) == Catch::Approx(std::max(x(0), 2.0)));
    CHECK(y(1) == Catch::Approx(std::max(x(1), 2.0)));
    CHECK(y(2) == Catch::Approx(std::max(x(2), 2.0)));
    CHECK(y(3) == Catch::Approx(std::max(x(3), 2.0)));
    CHECK(y(4) == Catch::Approx(std::max(x(4), 2.0)));
    CHECK(y(5) == Catch::Approx(std::max(x(5), 2.0)));

}

TEST_CASE("Elementwise min statically sized", "[utils]") {
    dlib::matrix<double, 6, 1> x;
    x =
            2.0, -1.0, 3.0, 5.0, 0.0, 1.5;

    dlib::matrix<double, 6, 1> y;
    y =
            cppsolnp::elementwise_min(x, 2.0);

    CHECK(y(0) == Catch::Approx(std::min(x(0), 2.0)));
    CHECK(y(1) == Catch::Approx(std::min(x(1), 2.0)));
    CHECK(y(2) == Catch::Approx(std::min(x(2), 2.0)));
    CHECK(y(3) == Catch::Approx(std::min(x(3), 2.0)));
    CHECK(y(4) == Catch::Approx(std::min(x(4), 2.0)));
    CHECK(y(5) == Catch::Approx(std::min(x(5), 2.0)));

}

TEST_CASE("Elementwise min dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(6);
    x =
            1.5, 2.0, -99.0, 15.0, 0.15, 1.5;

    dlib::matrix<double, 0, 1> y;
    y =
            cppsolnp::elementwise_min(x, 2.0);

    CHECK(y(0) == Catch::Approx(std::min(x(0), 2.0)));
    CHECK(y(1) == Catch::Approx(std::min(x(1), 2.0)));
    CHECK(y(2) == Catch::Approx(std::min(x(2), 2.0)));
    CHECK(y(3) == Catch::Approx(std::min(x(3), 2.0)));
    CHECK(y(4) == Catch::Approx(std::min(x(4), 2.0)));
    CHECK(y(5) == Catch::Approx(std::min(x(5), 2.0)));

}

TEST_CASE("Left vector min right vector max statically sized", "[utils]") {
    dlib::matrix<double, 3, 2> x;
    x =
            2.0, 1.0,
            5.0, 3.0,
            0.0, -1.5;

    dlib::matrix<double, 3, 2> y;
    y =
            cppsolnp::left_vector_min_right_vector_max(x);

    CHECK(y(0, 0) == Catch::Approx(std::min(x(0, 0), x(0, 1))));
    CHECK(y(0, 1) == Catch::Approx(std::max(x(0, 0), x(0, 1))));

    CHECK(y(1, 0) == Catch::Approx(std::min(x(1, 0), x(1, 1))));
    CHECK(y(1, 1) == Catch::Approx(std::max(x(1, 0), x(1, 1))));

    CHECK(y(2, 0) == Catch::Approx(std::min(x(2, 0), x(2, 1))));
    CHECK(y(2, 1) == Catch::Approx(std::max(x(2, 0), x(2, 1))));

}

TEST_CASE("Left vector min right vector max dynamically sized", "[utils]") {
    dlib::matrix<double> x(3, 2);
    x =
            2.0, 1.0,
            5.0, 3.0,
            0.0, -1.5;

    dlib::matrix<double> y(3, 2);
    y =
            cppsolnp::left_vector_min_right_vector_max(x);

    CHECK(y(0, 0) == Catch::Approx(std::min(x(0, 0), x(0, 1))));
    CHECK(y(0, 1) == Catch::Approx(std::max(x(0, 0), x(0, 1))));

    CHECK(y(1, 0) == Catch::Approx(std::min(x(1, 0), x(1, 1))));
    CHECK(y(1, 1) == Catch::Approx(std::max(x(1, 0), x(1, 1))));

    CHECK(y(2, 0) == Catch::Approx(std::min(x(2, 0), x(2, 1))));
    CHECK(y(2, 1) == Catch::Approx(std::max(x(2, 0), x(2, 1))));

}