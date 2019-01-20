

#include "stdafx.h"
#include <catch2/catch.hpp>

#include "utils.hpp"

TEST_CASE("Test equation solver", "[solve]") {
    dlib::matrix<double, 3, 3> A;
    A =
            2.0, 1.0, 3.0,
            2.0, 6.0, 8.0,
            6.0, 8.0, 18.0;

    dlib::matrix<double, 3, 1> y;

    y =
            1.0, 3.0, 5.0;

    dlib::matrix<double, 3, 1> x_sought;
    x_sought =
            3.0 / 10.0, 2.0 / 5.0, -0.0;

    dlib::matrix<double, 0, 1> x_answer;
    x_answer =
            cppsolnp::solve(A, y);

    REQUIRE(x_answer.nr() == x_sought.nr());

    CHECK(x_answer(0) == Approx(x_sought(0)));
    CHECK(x_answer(1) == Approx(x_sought(1)));
    CHECK(x_answer(2) == Approx(x_sought(2)));

}