#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <vector>
#include <functional>

#include <solver.hpp>


namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>);

template<typename RETURN_TYPE> using VectorFunction = std::function<RETURN_TYPE(std::vector<double>)>;

//template<typename M>
//VectorFunction<M> mapping_function(const cppsolnp::MatrixFunction<M> &f) {
//    return [f](dlib::matrix<double, 0, 1> param) {
//        std::vector<double> vector_param(param.begin(), param.end())
//        return f(vector_param);
//    };
//}

cppsolnp::MatrixFunction<double> objective_mapping_function(const VectorFunction<double> &f) {
    return [f](dlib::matrix<double, 0, 1> param) {
        std::vector<double> vector_param(param.begin(), param.end());
        return f(vector_param);
    };
}

cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> constraint_mapping_function(const VectorFunction<std::vector<double>> &f) {
    return [f](dlib::matrix<double, 0, 1> param) {
        std::vector<double> vector_param(param.begin(), param.end());
        std::vector<double> result = f(vector_param);
        return dlib::mat(result);
    };
}

double solve(
        std::vector<double> par_start_value,
        std::vector<double> par_upper_limit,
        std::vector<double> par_lower_limit,
        const std::function<double(std::vector<double>)> &obj_func,
        const std::function<std::vector<double>(std::vector<double>)> &eq_func,
        std::vector<double> eq_values,
        const std::function<std::vector<double>(std::vector<double>)> &ineq_func,
        std::vector<double> ineq_lower_bounds,
        std::vector<double> ineq_upper_bounds,
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



//    return 0.0;
    cppsolnp::SolverResult result = cppsolnp::solve_simple(
            dlib::mat(par_start_value),
            dlib::mat(par_upper_limit),
            dlib::mat(par_lower_limit),
            objective_mapping_function(obj_func),
            constraint_mapping_function(eq_func),
            dlib::mat(eq_values),
            constraint_mapping_function(ineq_func),
            dlib::mat(ineq_lower_bounds),
            dlib::mat(ineq_upper_bounds),
            rho,
            max_major_iter,
            max_minor_iter,
            delta,
            tolerance,
            debug);

    return result.solve_value;
}


//val_optim, x_optim, debug = pysolnp.solve(x0=start_point,
//                                          obj_func=alkyla_objective_function,
//                                          eq_func=alkyla_equality_function,
//                                          eq_values=equality_values,
//                                          ineq_func=alkyla_inequality_function,
//                                          ineq_lower_bounds=inequality_lower_bounds,
//                                          ineq_upper_bounds=inequality_upper_bounds,
//                                          rho=1.0,
//                                          max_major_iter=10,
//                                          max_minor_iter=10,
//                                          delta=1e-5,
//                                          tolerance=1e-4,
//                                          debug=True)


PYBIND11_MODULE(pysolnp, m
) {

m.doc() = R"pbdoc(

        -----------------------
        .. currentmodule:: pysolnp
        .. autosummary::
           :toctree: _generate
           solve
    )pbdoc";

m.def("solve", &solve, R"pbdoc(
        Solve using the SOLNP algorithm.
    )pbdoc",
    py::arg("par_start_value"),
    py::arg("par_upper_limit"),
    py::arg("par_lower_limit"),
    py::arg("obj_func"),
    py::arg("eq_func"),
    py::arg("eq_values"),
    py::arg("ineq_func"),
    py::arg("ineq_lower_bounds"),
    py::arg("ineq_upper_bounds"),
    py::arg("rho"),
    py::arg("max_major_iter"),
    py::arg("max_minor_iter"),
    py::arg("delta"),
    py::arg("tolerance"),
    py::arg("debug"));

    #ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
    #else
    m.attr("__version__") = "dev";
    #endif
}
