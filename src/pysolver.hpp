#ifndef CPP_SOLNP_PYSOLVER_HPP
#define CPP_SOLNP_PYSOLVER_HPP

#include <utility>
#include <vector>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "solnp.hpp"
#include "stdafx.h"
#include "utils.hpp"

#include <iostream>

namespace py = pybind11;

namespace pysolver {

    template<typename T>
    inline
    std::vector <T> py_list_to_std_vector(const py::list &list) {
        std::vector <T> result(list.size());
        int index = 0;
        for (auto &value : list) {
            result[index] = py::cast<T>(value);
            index++;
        }
        return result;
    }

    template<typename T>
    struct op_list_to_matrix {
        /*!
            This object defines a matrix expression that holds a std::vector.
            Thus it enables you to use a py::list as if it were a dlib::matrix.
        !*/
        explicit op_list_to_matrix(const py::list &list_) : list(py_list_to_std_vector<T>(list_)) {}

        const std::vector <T> list;

        // This expression wraps direct memory accesses so we use the lowest possible cost.
        const static long cost = 1;

        const static long NR = 0; // We don't know the length of the vector until runtime.  So we put 0 here.
        const static long NC = 1; // We do know that it only has one column (since it's a vector)
        typedef T type;
        // Since the py::list doesn't use a dlib memory manager we list the default one here.
        typedef dlib::default_memory_manager mem_manager_type;
        // The layout type also doesn't really matter in this case.  So we list row_major_layout
        // since it is a good default.
        typedef dlib::row_major_layout layout_type;

        typedef const T &const_ret_type;

        const_ret_type apply(long r, long) const {
            return list[r];
        }

        long nr() const { return list.size(); }

        long nc() const { return 1; }

        // This expression never aliases anything since it doesn't contain any matrix expression (it
        // contains only a std::vector which doesn't count since you can't assign a matrix expression
        // to a std::vector object).
        template<typename U>
        bool aliases(const dlib::matrix_exp <U> &) const { return false; }

        template<typename U>
        bool destructively_aliases(const dlib::matrix_exp <U> &) const { return false; }
    };

    template<typename T>
    dlib::matrix_op <op_list_to_matrix<T>> py_list_to_dlib_matrix(
            const py::list &list
    ) {
        typedef op_list_to_matrix<T> op;
        return dlib::matrix_op<op>(op(list));
    }

    template<typename T>
    inline
    py::list std_vector_to_py_list(
            const std::vector<T> &vector
    ) {
        py::list list(vector.size());
        for (auto row = 0L; row < vector.size(); ++row) {
            const T &value = vector(row);
            list[row] = value;
        }
        return list;
    }

    template<typename T>
    inline
    py::list dlib_matrix_to_py_list(
            const dlib::matrix<T, 0, 1> &matrix
    ) {
        py::list list(matrix.nr());
        for (auto row = 0L; row < matrix.nr(); ++row) {
            const T &value = matrix(row);
            list[row] = value;
        }
        return list;
    }

    struct Result {
        double solve_value;
        py::object optimum;
        int callbacks;

        Result(double solve_value, py::object optimum, int callbacks) :
                solve_value(solve_value), optimum(optimum), callbacks(callbacks) {
        }

    };


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
    );

}

namespace cppsolnp {
    template<typename RETURN_TYPE> using MatrixFunction = std::function<RETURN_TYPE(dlib::matrix<double, 0, 1>)>;

    struct SolverResult {
        double solve_value;
        dlib::matrix<double, 0, 1> optimum;
        cppsolnp::log_list_ptr log;
        int callbacks;

        SolverResult(double solve_value, dlib::matrix<double, 0, 1> optim, int function_calls,
                     cppsolnp::log_list_ptr log_list) :
                solve_value(solve_value), optimum(std::move(optim)), callbacks(function_calls), log(std::move(log_list)) {}
    };
}

#endif