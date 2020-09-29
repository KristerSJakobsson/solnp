
dlib::matrix<double, 4, 1> rosen_suzuki(const dlib::matrix<double, 4, 1> &m)
/*
Rosen-Suzuki function:
 This problem has only inequality constraints.
*/
{
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);
    const double x4 = m(3);


    // compute the rosen_suzuki function and return the result, equality constraint results and the inequality constraint results
    dlib::matrix<double, 4, 1> return_values(4);
    // Function value
    return_values(0) = x1 * x1 + x2 * x2 + 2 * x3 * x3 + x4 * x4 - 5 * x1 - 5 * x2 - 21 * x3 + 7 * x4;
    // Inequality constraints
    return_values(1) = 8 - x1 * x1 - x2 * x2 - x3 * x3 - x4 * x4 - x1 + x2 - x3 + x4;
    return_values(2) = 10 - x1 * x1 - 2 * x2 * x2 - x3 * x3 - 2 * x4 * x4 + x1 + x4;
    return_values(3) = 5 - 2 * x1 * x1 - x2 * x2 - x3 * x3 - 2 * x1 + x2 + x4;

    return return_values;
}


struct rosen_suzuki_functor {
public:
    rosen_suzuki_functor() = default;

    dlib::matrix<double, 4, 1> operator()(const dlib::matrix<double, 4, 1> &x) {
        return rosen_suzuki(x);
    }
};

