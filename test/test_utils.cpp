
#include <catch2/catch.hpp>

#include "../src/stdafx.h"
#include "../src/utils.hpp"


TEST_CASE("Conditional number statically sized", "[utils]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 1.0, 3.0,
            2.0, 6.0, 8.0,
            6.0, 8.0, 18.0;

    double cond =
            cppsolnp::conditional_number(A);

    CHECK(cond == Approx(31.413545501499051));

}

TEST_CASE("Conditional number dynamically sized", "[utils]") {
    dlib::matrix<double> A(3, 3);
    A =
            21.0, 14.0, -2.0,
            2.0, 1.0, -3.0,
            -21.0, 0.0, 18.0;

    double cond =
            cppsolnp::conditional_number(A);

    CHECK(cond == Approx(25.003922978899976));

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

    CHECK(x_divide_y(0) == Approx(x(0) / y(0)));
    CHECK(x_divide_y(1) == Approx(x(1) / y(1)));
    CHECK(x_divide_y(2) == Approx(x(2) / y(2)));

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

    CHECK(x_divide_y(0) == Approx(x(0) / y(0)));
    CHECK(x_divide_y(1) == Approx(x(1) / y(1)));
    CHECK(x_divide_y(2) == Approx(x(2) / y(2)));

}

TEST_CASE("Euclidean norm statically sized", "[utils]") {
    dlib::matrix<double, 4, 1> x;
    x =
            -2.0, 2.0, 3.0, 3.5;

    double norm = cppsolnp::euclidean_norm(x);

    CHECK(norm == Approx(5.408326913195984));

}

TEST_CASE("Euclidean norm dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(4);
    x =
            9.0, 15.0, 3.0, -2.5;

    double norm = cppsolnp::euclidean_norm(x);

    CHECK(norm == Approx(17.923448328934921));

}

TEST_CASE("Infinity norm statically sized", "[utils]") {
    dlib::matrix<double, 4, 1> x;
    x =
            -2.0, 2.0, 3.0, 3.5;

    double norm = cppsolnp::infinity_norm(x);

    CHECK(norm == Approx(3.5));

}

TEST_CASE("Infinity norm dynamically sized", "[utils]") {
    dlib::matrix<double, 0, 1> x(4);
    x =
            9.0, 15.0, 3.0, -2.5;

    double norm = cppsolnp::infinity_norm(x);

    CHECK(norm == Approx(15.0));

}
