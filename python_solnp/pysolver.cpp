#include "pysolver.hpp"


namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(pysolnp, m
) {
    py::class_<pysolver::Result>(m, "Result", py::dynamic_attr())
            .def(py::init<double, py::object, int, bool, py::object>())
            .def_readwrite("solve_value", &pysolver::Result::solve_value)
            .def_readwrite("optimum", &pysolver::Result::optimum)
            .def_readwrite("callbacks", &pysolver::Result::callbacks)
            .def_readwrite("converged", &pysolver::Result::converged)
            .def_readwrite("hessian_matrix", &pysolver::Result::hessian_matrix);

    m.doc() = R"pbdoc(

            -----------------------
            .. currentmodule:: pysolnp_test
            .. autosummary::
               :toctree: _generate
               solve
        )pbdoc";

    m.def("solve", &pysolver::python_solve, R"pbdoc(
            Solve using the SOLNP algorithm.
            )pbdoc",
          "obj_func"_a,
          "par_start_value"_a,
          "par_lower_limit"_a = py::cast<py::object>(Py_None),
          "par_upper_limit"_a = py::cast<py::object>(Py_None),
          "eq_func"_a = py::cast<py::object>(Py_None),
          "eq_values"_a = py::cast<py::object>(Py_None),
          "ineq_func"_a = py::cast<py::object>(Py_None),
          "ineq_lower_bounds"_a = py::cast<py::object>(Py_None),
          "ineq_upper_bounds"_a = py::cast<py::object>(Py_None),
          "rho"_a = 1.0,
          "max_major_iter"_a = 400,
          "max_minor_iter"_a = 800,
          "delta"_a = 1e-7,
          "tolerance"_a = 1e-8,
          "debug"_a = false);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}


namespace cppsolnp {

    CppsolnpResult solve_simple(const cppsolnp::MatrixFunction<double> &obj_func,
                              dlib::matrix<double, 0, 0> &parameter_data,
                              cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> eq_func,
                              std::shared_ptr<dlib::matrix<double, 0, 1>> eq_values,
                              cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> ineq_func,
                              std::shared_ptr<dlib::matrix<double, 0, 2>> ineq_data,
                              bool debug,
                              double rho, //penalty parameter
                              int maximum_major_iterations,
                              int maximum_minor_iterations,
                              const double &delta, // Step size in forward difference evaluation
                              const double &tolerance) {

        int function_calls = 0;

        cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> objective_function = [
                &obj_func,
                &eq_func,
                &eq_values,
                &ineq_func,
                &ineq_data,
                &function_calls](const dlib::matrix<double, 0, 1> &point) {

            auto number_of_parameters = 1L;

            if (eq_func && eq_values) {
                number_of_parameters += (*eq_values).nr();
            }

            if (ineq_func && ineq_data) {
                number_of_parameters += (*ineq_data).nr();
            }

            dlib::matrix<double, 0, 1> result = dlib::matrix<double, 0, 1>(number_of_parameters);

            double optimal_function_value = obj_func(point);
            result(0) = optimal_function_value;

            auto save_to_index = 1L;

            if (eq_func && eq_values) {
                dlib::matrix<double, 0, 1> equality_function_value = eq_func(point);

                if (equality_function_value.nr() != (*eq_values).nr()) {
                    throw std::invalid_argument(
                            "Equality function evaluated to a different length than the equality values.");
                }

                // Equality constraints
                for (auto row = 0L; row < (*eq_values).nr(); row++) {
                    result(save_to_index) =
                            equality_function_value(row) - (*eq_values)(row); // Offset by right-hand value to get == 0
                    save_to_index++;
                }

            }

            if (ineq_func && ineq_data) {
                dlib::matrix<double, 0, 1> inequality_function_value = ineq_func(point);

                if (inequality_function_value.nr() != (*ineq_data).nr()) {
                    throw std::invalid_argument(
                            "Inequality function evaluated to a different length than the inequality bounds.");
                }

                // Inequality constraints
                for (auto row = 0L; row < (*ineq_data).nr(); row++) {
                    result(save_to_index) = inequality_function_value(row);
                    save_to_index++;
                }

            }

            function_calls += 1;
            return result;
        };

        std::shared_ptr<std::vector<std::string>> logger;
        if (debug) {
            logger = std::make_shared<std::vector<std::string>>(std::vector<std::string>());
        }

        /*
         * No Hessian Matrix provided means it assumes the unit matrix
         */

        dlib::matrix<double, 0, 0> inequality_limits;
        if (ineq_data) {
            inequality_limits = *ineq_data;
        }

        cppsolnp::SolveResult result = cppsolnp::solnp(
                objective_function,
                parameter_data,
                inequality_limits,
                logger,
                rho,
                maximum_major_iterations,
                maximum_minor_iterations,
                delta,
                tolerance);

        CppsolnpResult final_result(result.solve_value, result.optimum, function_calls, result.converged, logger, result.hessian_matrix);

        return final_result;
    }
}

namespace pysolver {


    cppsolnp::MatrixFunction<double> objective_mapping_function(const py::function &f) {
        auto python_function = f.cast<py::function>();
        return [python_function](const dlib::matrix<double, 0, 1> &param) {
            const py::list &list_param = dlib_1d_matrix_to_py_list<double>(param);
            double result = py::float_(python_function(list_param));
            return result;
        };
    }

    cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>>
    constraint_mapping_function(const py::function &f, const dlib::matrix<double> &val, const std::string &name) {
        auto python_function = f.cast<py::function>();
        return [python_function, val, name](const dlib::matrix<double, 0, 1> &param) {
            const py::list &list_param = dlib_1d_matrix_to_py_list<double>(param);
            const py::list &result = python_function(list_param);
            if (result.size() != (size_t) val.nr()) {
                std::string expected_size = std::to_string(val.nr());
                std::string actual_size = std::to_string(result.size());
                throw std::invalid_argument(
                        "The " + name + " callback function returned a value of length " + actual_size +
                        " but comparison values are of length " + expected_size);
            }
            const dlib::matrix<double, 0, 1> &result_matrix = py_list_to_dlib_matrix<double>(result);
            return result_matrix;
        };
    }

    Result python_solve(
            const py::function &obj_func,
            const py::list &par_start_value,
            const py::object &par_lower_limit,
            const py::object &par_upper_limit,
            const py::object &eq_func,
            const py::object &eq_values,
            const py::object &ineq_func,
            const py::object &ineq_lower_bounds,
            const py::object &ineq_upper_bounds,
            double rho,
            int max_major_iter,
            int max_minor_iter,
            double delta,
            double tolerance,
            bool debug
    ) {

        dlib::matrix<double, 0, 0> parameter_data;
        if (!par_lower_limit.is_none() && !par_upper_limit.is_none()) {
            parameter_data = dlib::join_rows(pysolver::py_list_to_dlib_matrix<double>(par_start_value),
                                             dlib::join_rows(pysolver::py_list_to_dlib_matrix<double>(par_lower_limit),
                                                             pysolver::py_list_to_dlib_matrix<double>(
                                                                     par_upper_limit)));
        } else if (par_lower_limit.is_none() && par_upper_limit.is_none()) {
            parameter_data = pysolver::py_list_to_dlib_matrix<double>(par_start_value);
        } else {
            throw std::invalid_argument(
                    "Bad input: Can not only provide equality or inequaltiy constraints, please provide both.");
        }

        cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> equality_function(nullptr);
        std::shared_ptr<dlib::matrix<double, 0, 1>> equality_function_values(nullptr);
        if (!eq_func.is_none() && !eq_values.is_none()) {
            equality_function_values = std::make_shared<dlib::matrix<double, 0, 1
            >>(pysolver::py_list_to_dlib_matrix<double>(eq_values));
            equality_function = constraint_mapping_function(eq_func, *equality_function_values, "equality");
        } else if (eq_func.is_none() && eq_values.is_none()) {
            // Do nothing
        } else {
            throw std::invalid_argument(
                    "Bad input: Must provide equality function together with equality values or not provide one at all.");
        }

        cppsolnp::MatrixFunction<dlib::matrix<double, 0, 1>> inequality_function(nullptr);
        std::shared_ptr<dlib::matrix<double, 0, 2>> inequality_function_data(nullptr);
        if (!ineq_func.is_none() && !ineq_lower_bounds.is_none() && !ineq_upper_bounds.is_none()) {
            inequality_function_data = std::make_shared<dlib::matrix<double, 0, 2 >>(
                    dlib::join_rows(pysolver::py_list_to_dlib_matrix<double>(ineq_lower_bounds),
                                    pysolver::py_list_to_dlib_matrix<double>(ineq_upper_bounds)));
            inequality_function = constraint_mapping_function(ineq_func, *inequality_function_data, "inequality");
        } else if (ineq_func.is_none() && ineq_lower_bounds.is_none() && ineq_upper_bounds.is_none()) {
            // Do nothing
        } else {
            throw std::invalid_argument(
                    "Bad input: Must provide inequality function together with upper and lower bounds or not provide one at all.");
        }

        cppsolnp::CppsolnpResult result = cppsolnp::solve_simple(
                objective_mapping_function(obj_func),
                parameter_data,
                equality_function,
                equality_function_values,
                inequality_function,
                inequality_function_data,
                debug,
                rho,
                max_major_iter,
                max_minor_iter,
                delta,
                tolerance);

        if (result.log) {
            for (const auto& val : *result.log) {
                std::cout << val << std::endl;
            }
        }

        const py::object &return_optimum = pysolver::dlib_1d_matrix_to_py_list<double>(result.optimum);
        const py::object &return_hessian_matrix = pysolver::dlib_2d_matrix_to_py_nested_list<double>(result.hessian_matrix);

        pysolver::Result return_value(result.solve_value, return_optimum, result.callbacks, result.converged, return_hessian_matrix);
        return return_value;
    }

}


