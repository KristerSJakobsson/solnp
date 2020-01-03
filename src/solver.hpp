#ifndef CPP_SOLNP_SOLVER_HPP
#define CPP_SOLNP_SOLVER_HPP

#include <utility>

#include "solnp.hpp"

/*
 * This file holds the entry functions for using CPP SOLNP
 */


namespace cppsolnp {

    template<typename RETURN_TYPE> using MatrixFunction = std::function<RETURN_TYPE(dlib::matrix<double, 0, 1>)>;

    struct SolverResult {
        double solve_value;
        dlib::matrix<double, 0, 1> optimum;
        cppsolnp::log_list log;
        int callbacks;

        SolverResult(double solve_value, dlib::matrix<double, 0, 1> optim, int function_calls, const log_list_ptr &log_list) :
                solve_value(solve_value), optimum(std::move(optim)), callbacks(function_calls), log(*log_list)  {}
    };

    SolverResult solve_simple(const dlib::matrix<double, 0, 1> &par_start,
                              const dlib::matrix<double, 0, 1> &par_lower_bound,
                              const dlib::matrix<double, 0, 1> &par_upper_bound,
                              MatrixFunction<double> obj_func,
                              MatrixFunction<dlib::matrix<double, 0, 1>> eq_func,
                              const dlib::matrix<double, 0, 1> &eq_values,
                              MatrixFunction<dlib::matrix<double, 0, 1>> ineq_func,
                              const dlib::matrix<double, 0, 1> &ineq_lower_bounds,
                              const dlib::matrix<double, 0, 1> &ineq_upper_bounds,
                              bool debug, // TODO: Wire this up
                              double rho = 1.0, //penalty parameter
                              int maximum_major_iterations = 10,
                              int maximum_minor_iterations = 10,
                              const double &delta = 1e-5, // Step size in forward difference evaluation
                              const double &tolerance = 1e-4) {

        if (par_start.nr() != par_lower_bound.nr() ||
            par_start.nc() != par_upper_bound.nc()) {
            throw std::invalid_argument("Bad input: Length of parameter lower bound must match length of upper bound.");
        }

        dlib::matrix<double, 0, 3> parameter_data = dlib::join_rows(par_start,
                                                                    dlib::join_rows(par_lower_bound, par_upper_bound));

        if (ineq_lower_bounds.nr() != ineq_upper_bounds.nr() ||
            ineq_lower_bounds.nc() != ineq_upper_bounds.nc()) {
            throw std::invalid_argument(
                    "Bad input: Length of inequality constraint lower bound must match length of upper bound.");
        }

        dlib::matrix<double, 0, 2> inequality_constraint_data = dlib::join_rows(ineq_lower_bounds, ineq_upper_bounds);

        int function_calls = 0;

        MatrixFunction<dlib::matrix<double, 0, 1>> objective_function = [
                &obj_func,
                &eq_func,
                &eq_values,
                &ineq_func,
                &inequality_constraint_data,
                &function_calls](const dlib::matrix<double, 0, 1> &point) {

            dlib::matrix<double, 0, 1> result = dlib::matrix<double, 0, 1>(
                    1 + eq_values.nr() + inequality_constraint_data.nr());

            double optimal_function_value = obj_func(point);
            dlib::matrix<double, 0, 1> equality_function_value = eq_func(point);
            dlib::matrix<double, 0, 1> inequality_function_value = ineq_func(point);

            if (equality_function_value.nr() != eq_values.nr()) {
                throw std::invalid_argument(
                        "Error: Equality function evaluated to a different length than the equality values.");
            }

            if (inequality_function_value.nr() != inequality_constraint_data.nr()) {
                throw std::invalid_argument(
                        "Error: Inequality function evaluated to a different length than the inequality bounds.");
            }

            result(0) = optimal_function_value;
            // Equality constraints
            for (auto row = 0L; row < eq_values.nr(); row++) {
                result(1 + row) =
                        equality_function_value(row) - eq_values(row); // Offset by right-hand value to get == 0
            }
            // Inequality constraints
            for (auto row = 0L; row < inequality_constraint_data.nr(); row++) {
                result(1 + eq_values.nr() + row) = inequality_function_value(row);
            }
            function_calls += 1;
            return result;
        };

        cppsolnp::log_list_ptr logger(new cppsolnp::log_list());

        /*
         * No Hessian Matrix provided means it assumes the unit matrix
         */

        double result = cppsolnp::solnp(
                objective_function,
                parameter_data,
                inequality_constraint_data,
                logger,
                nullptr,
                rho,
                maximum_major_iterations,
                maximum_minor_iterations,
                delta,
                tolerance);

        dlib::matrix<double, 0, 1> final_vector = dlib::colm(parameter_data, 0);

        SolverResult final_result(result, final_vector, function_calls, logger);

        return final_result;
    }


    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    double solve(
            functor_model functor,
            parameter_input &parameter_data,
            const inequality_constraint_vectors &inequality_constraint_data,
            const cppsolnp::log_list_ptr &event_log = nullptr,
            // Below are control variables.
            double rho = 1.0, //penalty parameter
            int maximum_major_iterations = 10,
            int maximum_minor_iterations = 10,
            const double &delta = 1e-5, // Step size in forward difference evaluation
            const double &tolerance = 1e-4
    ) {
        /*
         * No Hessian Matrix provided means it assumes the unit matrix
         */
        return cppsolnp::solnp(
                functor,
                parameter_data,
                inequality_constraint_data,
                event_log,
                nullptr,
                rho,
                maximum_major_iterations,
                maximum_minor_iterations,
                delta,
                tolerance);
    }


    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    double solve(
            functor_model functor,
            parameter_input &parameter_data,
            const inequality_constraint_vectors &inequality_constraint_data,
            const dlib::matrix<double> hessian_matrix,
            const cppsolnp::log_list_ptr &event_log = nullptr,
            // Below are control variables.
            double rho = 1.0, //penalty parameter
            int maximum_major_iterations = 10,
            int maximum_minor_iterations = 10,
            const double &delta = 1e-5, // Step size in forward difference evaluation
            const double &tolerance = 1e-4
    ) {
        /*
         * No Hessian Matrix provided -> pass nullptr to solnp
         */
        return cppsolnp::solnp(
                functor,
                parameter_data,
                inequality_constraint_data,
                event_log,
                &hessian_matrix,
                rho,
                maximum_major_iterations,
                maximum_minor_iterations,
                delta,
                tolerance);
    }
} // namespace cppsolnp

#endif // CPP_SOLNP_SOLVER_HPP