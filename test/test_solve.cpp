

#include "../src/stdafx.h"
#include <catch2/catch.hpp>

#include "../src/solve.hpp"

TEST_CASE("Determined system", "[solve]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 1.0, 3.0,
            2.0, 6.0, 8.0,
            6.0, 8.0, 18.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    dlib::matrix<double, 0, 1> y_calculated = A * x_answer;

    CHECK(y_calculated(0) == Approx(y(0)));
    CHECK(y_calculated(1) == Approx(y(1)));
    CHECK(y_calculated(2) == Approx(y(2)));

}

TEST_CASE("Undetermined system", "[solve]") {
    dlib::matrix<double, 2, 3> A;
    A =
            1.0, 2.0, 0.0,
            0.0, 4.0, 3.0;

    dlib::matrix<double, 2, 1> y;

    y =
            8.0, 18.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    dlib::matrix<double, 0, 1> y_calculated = A * x_answer;

    CHECK(y_calculated(0) == Approx(y(0)));
    CHECK(y_calculated(1) == Approx(y(1)));

}


TEST_CASE("Upper triangular matrix", "[solve]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 1.0, 3.0,
            0.0, 6.0, 8.0,
            0.0, 0.0, 18.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    dlib::matrix<double, 0, 1> y_calculated = A * x_answer;

    CHECK(y_calculated(0) == Approx(y(0)));
    CHECK(y_calculated(1) == Approx(y(1)));
    CHECK(y_calculated(2) == Approx(y(2)));

}


TEST_CASE("Lower triangular matrix", "[solve]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 0.0, 0.0,
            2.0, 6.0, 0.0,
            6.0, 8.0, 18.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    dlib::matrix<double, 0, 1> y_calculated = A * x_answer;

    CHECK(y_calculated(0) == Approx(y(0)));
    CHECK(y_calculated(1) == Approx(y(1)));
    CHECK(y_calculated(2) == Approx(y(2)));

}


TEST_CASE("Diagonal matrix", "[solve]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 0.0, 0.0,
            0.0, 6.0, 0.0,
            0.0, 0.0, 18.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    dlib::matrix<double, 0, 1> y_calculated = A * x_answer;

    CHECK(y_calculated(0) == Approx(y(0)));
    CHECK(y_calculated(1) == Approx(y(1)));
    CHECK(y_calculated(2) == Approx(y(2)));

}