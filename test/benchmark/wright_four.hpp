
dlib::matrix<double, 4, 1> wright_four(const dlib::matrix<double, 5, 1> &m)
/*
Wright4 function from the original documentation of SOLNP:
 Wright4: this problem has only equality constraints
*/
{
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);
    const double x4 = m(3);
    const double x5 = m(4);


    // compute the wright_four function and return the result, equality constraint results and teh inequality constraint results
    dlib::matrix<double, 4, 1> return_values(4);
    // Function value
    return_values(0) = (x1 - 1) * (x1 - 1) +
                       (x1 - x2) * (x1 - x2) +
                       (x2 - x3) * (x2 - x3) * (x2 - x3) +
                       (x3 - x4) * (x3 - x4) * (x3 - x4) * (x3 - x4) +
                       (x4 - x5) * (x4 - x5) * (x4 - x5) * (x4 - x5);
    // Equality constraints
    return_values(1) = x1 + x2 * x2 + x3 * x3 * x3 - 2 - 3 * std::sqrt(2.0);
    return_values(2) = x2 - x3 * x3 + x4 + 2 - 2 * std::sqrt(2);
    return_values(3) = x1 * x5 - 2;
    return return_values;
}


struct wright_four_functor {
public:
    wright_four_functor() = default;;

    dlib::matrix<double, 4, 1> operator()(const dlib::matrix<double, 5, 1> &x) {
        return wright_four(x);
    }
};

