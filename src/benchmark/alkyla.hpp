#ifndef CPP_SOLNP_BENCHMARK_ALKYLA_HPP
#define CPP_SOLNP_BENCHMARK_ALKYLA_HPP

dlib::matrix<double, 8, 1> alkyla(const dlib::matrix<double, 10, 1> &m)
/*
Alkyla function from the original documentation of SOLNP:
 ALKYLA: this problem has both equality and inequality constraints, as well as the parameter bounds.
*/
{
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


    // compute the alkyla function and return the result, equality constraint results and the inequality constraint results
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
    alkyla_functor() = default;;

    dlib::matrix<double, 8, 1> operator()(const dlib::matrix<double, 10, 1> &x) {
        return alkyla(x);
    }
};

#endif