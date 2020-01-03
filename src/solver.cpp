#include <utility>
#include <vector>
#include <functional>

#include "solver.hpp"

template<typename RETURN_TYPE> using VectorFunction = std::function<RETURN_TYPE(std::vector<double>)>;

// Map between std::vector <-> dlib::matrix
cppsolnp::MatrixFunction<double> objective_mapping_function(const VectorFunction<double> &f) {
    return [f](dlib::matrix<double, 0, 1> param) {
        std::vector<double> vector_param(param.begin(), param.end());
        double result = f(vector_param);
        return result;
    };
}

cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>>
constraint_mapping_function(const VectorFunction<std::vector<double>> &f) {
    return [f](dlib::matrix<double, 0, 1> param) {
        std::vector<double> vector_param(param.begin(), param.end());
        dlib::matrix<double, 0, 1> result = dlib::mat(f(vector_param));
        return result;
    };

}

struct Result {
    double solve_value;
    std::vector<double> optimum;
    int callbacks;
    cppsolnp::log_list log;

    Result(double solve_value, std::vector<double> optimum, cppsolnp::log_list log, int callbacks) :
            solve_value(solve_value), optimum(std::move(optimum)), log(std::move(log)), callbacks(callbacks) {}
};

Result solve(
        std::vector<double> &par_start_value,
        std::vector<double> &par_lower_limit,
        std::vector<double> &par_upper_limit,
        const std::function<double(std::vector<double>)> &obj_func,
        const std::function<std::vector<double>(std::vector<double>)> &eq_func,
        std::vector<double> &eq_values,
        const std::function<std::vector<double>(std::vector<double>)> &ineq_func,
        std::vector<double> &ineq_lower_bounds,
        std::vector<double> &ineq_upper_bounds,
        double rho,
        int max_major_iter,
        int max_minor_iter,
        double delta,
        double tolerance,
        bool debug
) {


//    VectorFunction<double> mapping_function(const cppsolnp::MatrixFunction<double> &f) {
//        return [f](dlib::matrix<double, 0, 1> param) {
//            std::vector<double> vector_param(param.begin(), param.end())
//            return f(vector_param);
//        };
//    }
//
//    std::cout << par_start_value << std::endl;
//    std::cout << par_lower_limit << std::endl;
//    std::cout << par_upper_limit << std::endl;
//    std::cout << eq_values << std::endl;
//    std::cout << ineq_lower_bounds << std::endl;
//    std::cout << ineq_upper_bounds << std::endl;

//    return 0.0;
    cppsolnp::SolverResult result = cppsolnp::solve_simple(
            dlib::mat(par_start_value),
            dlib::mat(par_lower_limit),
            dlib::mat(par_upper_limit),
            objective_mapping_function(obj_func),
            constraint_mapping_function(eq_func),
            dlib::mat(eq_values),
            constraint_mapping_function(ineq_func),
            dlib::mat(ineq_lower_bounds),
            dlib::mat(ineq_upper_bounds),
            debug,
            rho,
            max_major_iter,
            max_minor_iter,
            delta,
            tolerance);

    std::vector<double> return_optimum(result.optimum.begin(), result.optimum.end());
    Result return_value(result.solve_value, return_optimum, result.log, result.callbacks);
    return return_value;
}
