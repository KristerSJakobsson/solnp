
dlib::matrix<double, 2, 1> box(const dlib::matrix<double, 3, 1> &m)
/*
The Box function from the original documentation of SOLNP:
 BOX: this problem has one equality constraints and variable bounds
*/
{
    const double x1 = m(0);
    const double x2 = m(1);
    const double x3 = m(2);


    // compute the alkyla function and return the result, equality constraint results and teh inequality constraint results
    dlib::matrix<double, 2, 1> return_values(2);
    // Function value
    return_values(0) = -1 * x1 * x2 * x3;
    // Equality constraints
    return_values(1) = 4 * x1 * x2 + 2 * x2 * x3 + 2 * x3 * x1 - 100;
    return return_values;
}


struct box_functor {
public:
    box_functor() = default;;

    dlib::matrix<double, 2, 1> operator()(const dlib::matrix<double, 3, 1> &x) {
        return box(x);
    }
};

