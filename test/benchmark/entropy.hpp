
double inline vector_entropy(const dlib::matrix<double, 0, 1> &vector) {
    double result = 0.0;
    for (auto value : vector) {
        result -= std::log(value);
    }
    dlib::matrix<double, 0, 0> ones_matrix = dlib::ones_matrix<double>(vector.nr(), vector.nc());
    result -= std::log(cppsolnp::euclidean_norm(vector - ones_matrix) + 0.1);
    return result;
}

dlib::matrix<double, 0, 1> entropy(const dlib::matrix<double, 0, 1> &m)
/*
Entropy function from the original documentation of SOLNP:
 ENTROPY: this is a nonconvex problem
*/
{
    // compute the entropy function and return the result, equality constraint results and teh inequality constraint results
    dlib::matrix<double, 2, 1> return_values(2);
    // Function value
    return_values(0) = vector_entropy(m);
    // Equality constraints
    return_values(1) = dlib::sum(m) - m.nr();
    return return_values;
}


struct entropy_functor {
public:
    entropy_functor() = default;

    dlib::matrix<double, 2, 1> operator()(const dlib::matrix<double, 0, 1> &x) {
        return entropy(x);
    }
};

