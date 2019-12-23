
dlib::matrix<double, 4, 1> wright_nine(const dlib::matrix<double, 5, 1> &m)
/*
Wright9 function from the original documentation of SOLNP:
 Wright9: this problem has inequality constraints.
*/
{
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);
    const double x4 = m(3);
    const double x5 = m(4);


    // compute the wright_nine function and return the result, equality constraint results and teh inequality constraint results
    dlib::matrix<double, 4, 1> return_values(4);
    // Function value
    return_values(0) = 10 * x1 * x4 - 6 * x3 * x2 * x2 + x2 * x1 * x1 * x1 + 9 * std::sin(x5 - x3) +
                       x5 * x5 * x5 * x5 * x4 * x4 * x2 * x2 * x2;
    // Inequality constraints
    return_values(1) = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5;
    return_values(2) = x3 * x1 * x1 + x4 * x5;
    return_values(3) = x4 * x2 * x2 + 10 * x1 * x5;
    return return_values;
}


struct wright_nine_functor {
public:
    wright_nine_functor() = default;;

    dlib::matrix<double, 4, 1> operator()(const dlib::matrix<double, 5, 1> &x) {
        return wright_nine(x);
    }
};

