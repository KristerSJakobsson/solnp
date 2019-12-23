
dlib::matrix<double, 4, 1> powell(const dlib::matrix<double, 5, 1> &m)
/*
This function computes what the Powell function from the
original documentaiton of SOLNP.
*/
{
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
    powell_functor() = default;;

    dlib::matrix<double, 4, 1> operator()(const dlib::matrix<double, 5, 1> &x) {
        return powell(x);
    }
};
