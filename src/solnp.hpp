#ifndef CPP_SOLNP_SOLNP_HPP
#define CPP_SOLNP_SOLNP_HPP


#include "stdafx.h"

#include "utils.hpp"
#include "subnp.hpp"

/* This is an interior method QP solver. */
namespace cppsolnp {

    template<typename T> bool isfinite(T arg)
    {
        return arg == arg &&
               arg != std::numeric_limits<T>::infinity() &&
               arg != -std::numeric_limits<T>::infinity();
    }

    struct SolveResult {
        double solve_value;
        dlib::matrix<double, 0, 1> optimum;
        bool converged;
        dlib::matrix<double> hessian_matrix;

        SolveResult(double solve_value, dlib::matrix<double, 0, 1> optimum, bool converged, dlib::matrix<double> hessian_matrix) :
                solve_value(solve_value), optimum(optimum), converged(converged), hessian_matrix(hessian_matrix) {}
    };

    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    SolveResult solnp(
            functor_model functor,
            parameter_input &parameter_data,
            const inequality_constraint_vectors &inequality_constraint_data,
            dlib::matrix<double> &hessian_matrix,
            // Optional input
            const std::shared_ptr<std::vector<std::string>> &event_log = nullptr,
            double rho = 1.0, //penalty parameter
            int maximum_major_iterations = 400,
            int maximum_minor_iterations = 800,
            const double &delta = 1e-7, // Step size in forward difference evaluation
            const double &tolerance = 1e-8
    ) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<parameter_input>::value);
        COMPILE_TIME_ASSERT(dlib::is_matrix<inequality_constraint_vectors>::value);
        COMPILE_TIME_ASSERT(parameter_input::NC <= 3 && parameter_input::NC >= 0);
        COMPILE_TIME_ASSERT(inequality_constraint_vectors::NC <= 3);

        if (!isfinite<double>(tolerance * tolerance)) {
            throw std::invalid_argument("Tolerance set too low.");
        }

        /* Split the parameter from the
        parameter constraints if applicable.*/
        dlib::matrix<double> parameters;
        dlib::matrix<double> parameter_bounds;
        dlib::matrix<double> objective_function_gradient;

        auto number_of_parameters = parameter_data.nr(),
                parameter_vector_width = parameter_data.nc();

        std::pair<bool, bool> lagrangian_parameters_bounded;
        lagrangian_parameters_bounded.first = true;
        lagrangian_parameters_bounded.second = false;

        if (parameter_vector_width == 1) {
            parameters = dlib::colm(parameter_data, 0);
            lagrangian_parameters_bounded.first = false;
        } else if (parameter_vector_width == 2) {
            parameters = 0.5 * (dlib::colm(parameter_data, 0) + dlib::colm(parameter_data, 1));
            parameter_bounds = dlib::colm(parameter_data, dlib::range(0, 1));
        } else if (parameter_vector_width == 3) {
            parameters = dlib::colm(parameter_data, 0);
            parameter_bounds = dlib::colm(parameter_data, dlib::range(1, 2));

        } else {
            throw std::invalid_argument("Parameter array must have three columns or less.");
        }

        if (lagrangian_parameters_bounded.first) {
            if (parameter_vector_width == 2)
            {
                if (dlib::min(dlib::colm(parameter_data, 1) - dlib::colm(parameter_data, 0)) <= 0) {
                    throw std::invalid_argument(
                            "The lower bounds of the parameter constraints must be strictly less than the upper bounds.");
                }
            }
            else if (parameter_vector_width == 3)
            {
                if (dlib::min(dlib::colm(parameter_data, 2) - dlib::colm(parameter_data, 1)) <= 0) {
                    throw std::invalid_argument(
                            "The lower bounds of the parameter constraints must be strictly less than the upper bounds.");
                }
                if (dlib::min(parameters - dlib::colm(parameter_data, 1)) <= 0 ||
                    dlib::min(dlib::colm(parameter_data, 2) - parameters) <= 0)
                {
                    throw std::invalid_argument("Initial parameter values must be within the bounds.");
                }
            }
        }
        /* Assume inequality constraints exist for now.
        TODO: Either properly handle it or just use specialization.
        Lots of possibilities, little time.*/


        int inequality_constraints_vector_length = inequality_constraint_data.nr(),
                inequality_constraints_vector_width = inequality_constraint_data.nc();

        int number_of_inequality_constraints;

        dlib::matrix<double, inequality_constraint_vectors::NR, 1> temporary_inequality_guess; // ib0
        dlib::matrix<double, inequality_constraint_vectors::NR, 2> temporary_inequality_constraints; //ib


        if (inequality_constraints_vector_width == 0) {
            number_of_inequality_constraints = 0;
        } else {
            number_of_inequality_constraints = inequality_constraints_vector_length;


            if (inequality_constraints_vector_width == 3) {
                temporary_inequality_guess = dlib::colm(inequality_constraint_data, 0); // ib0
                temporary_inequality_constraints = dlib::colm(inequality_constraint_data, dlib::range(1, 2)); // ib

                if (dlib::min(temporary_inequality_guess - dlib::colm(temporary_inequality_constraints, 0)) <= 0 ||
                    dlib::min(dlib::colm(temporary_inequality_constraints, 1) - temporary_inequality_guess) <= 0) {
                    throw std::invalid_argument("Initial inequalities must be within bounds.");
                }

            } else if (inequality_constraints_vector_width == 2) {
                if (dlib::min(dlib::colm(inequality_constraint_data, 1) - dlib::colm(inequality_constraint_data, 0)) <=
                    0) {
                    throw std::invalid_argument(
                            "The lower bounds of the inequality constraints must be strictly less than the upper bounds.");
                }
                temporary_inequality_guess =
                        0.5 * (dlib::colm(inequality_constraint_data, 0) + dlib::colm(inequality_constraint_data, 1));
                temporary_inequality_constraints = inequality_constraint_data;
            } else if (inequality_constraints_vector_width == 1) {
                number_of_inequality_constraints = 0;
            } else {
                throw std::invalid_argument("Inequality constraints must have 2 or 3 columns.");
            }
            if (number_of_inequality_constraints > 0) {
                if (lagrangian_parameters_bounded.first) {
                    parameter_bounds = dlib::join_cols(temporary_inequality_constraints, parameter_bounds);
                } else {
                    parameter_bounds = temporary_inequality_constraints;
                }
                parameters = dlib::join_cols(temporary_inequality_guess, parameters);
            }
        }

        // Here we could release the temporary matrixes.
        if (hessian_matrix.nr() != number_of_parameters + number_of_inequality_constraints ||
            hessian_matrix.nc() != number_of_parameters + number_of_inequality_constraints) {
            throw std::invalid_argument("The provided hessian matrix override was of invalid dimension.");
        }

        if (lagrangian_parameters_bounded.first || number_of_inequality_constraints > 0) {
            lagrangian_parameters_bounded.second = true;
        }

        dlib::matrix<double> cost_vector;
        cost_vector = functor(dlib::rowm(parameters, dlib::range(number_of_inequality_constraints,
                                                                 number_of_inequality_constraints +
                                                                 number_of_parameters - 1))); // ob

        if (event_log)
            event_log->push_back("Updated parameters: " +
                                 to_string(dlib::rowm(parameters, dlib::range(number_of_inequality_constraints,
                                                                              number_of_inequality_constraints +
                                                                              number_of_parameters - 1)),
                                           true));

        auto cost_vector_length = cost_vector.nr(),
                cost_vector_width = cost_vector.nc();
        if (cost_vector_width != 1) {
            throw std::invalid_argument("sqp_min cost function must return only 1 column.");
        }
        if (cost_vector_length < number_of_inequality_constraints + 1) {
            throw std::invalid_argument(
                    "sqp_min the number of constraints in the cost function does not match the call to sqp_min.");
        }

        auto number_of_equality_constraints = cost_vector_length - 1 - number_of_inequality_constraints;
        auto number_of_constraints = cost_vector_length - 1;

        double objective_function_value = cost_vector(0);
        std::vector<double> objective_function_value_history;
        objective_function_value_history.push_back(objective_function_value);
        dlib::matrix<double, 3, 1> t;
        t = dlib::zeros_matrix<double>(3, 1);

        dlib::matrix<double, 0, 1> lagrangian_multipliers;
        dlib::matrix<double, 0, 1> constraints;

        if (number_of_constraints != 0) {
            lagrangian_multipliers = dlib::zeros_matrix<double>(number_of_constraints, 1);

            constraints = dlib::rowm(cost_vector, dlib::range(1, number_of_constraints));

            if (number_of_inequality_constraints != 0) {
                if (dlib::min(
                        dlib::rowm(constraints,
                                   dlib::range(number_of_equality_constraints, number_of_constraints - 1)) -
                        dlib::subm(parameter_bounds, dlib::rectangle(0, 0, 0, number_of_inequality_constraints - 1))
                ) > 0.0 &&
                    dlib::min(
                            dlib::subm(parameter_bounds,
                                       dlib::rectangle(1, 0, 1, number_of_inequality_constraints - 1)) -
                            dlib::rowm(constraints,
                                       dlib::range(number_of_equality_constraints, number_of_constraints - 1))
                    ) > 0.0) {
                    dlib::set_rowm(parameters, dlib::range(0, number_of_inequality_constraints - 1)) =
                            dlib::rowm(constraints,
                                       dlib::range(number_of_equality_constraints, number_of_constraints - 1));
                }
                dlib::set_rowm(constraints, dlib::range(number_of_equality_constraints, number_of_constraints - 1)) =
                        dlib::rowm(constraints,
                                   dlib::range(number_of_equality_constraints, number_of_constraints - 1)) -
                        dlib::rowm(parameters, dlib::range(0, number_of_inequality_constraints - 1));
            }

            t(1) = euclidean_norm(constraints);

            if (std::max<double>(t(1) - 10.0 * tolerance, number_of_inequality_constraints) <= 0) {
                rho = 0.0;
            }

        } else {
            lagrangian_multipliers = dlib::zeros_matrix<double>(1, 1); // lambda = 0
        }

        double mu = number_of_parameters;
        int iteration = 0;

        subnp<functor_model> sub_problem(
                functor,
                number_of_parameters,
                number_of_equality_constraints,
                number_of_inequality_constraints,
                lagrangian_parameters_bounded);

        while (iteration < maximum_major_iterations) {
            ++iteration;

            sub_problem(parameters,
                        parameter_bounds,
                        lagrangian_multipliers,
                        cost_vector,
                        hessian_matrix,
                        mu,
                        rho,
                        maximum_minor_iterations,
                        delta,
                        tolerance,
                        event_log);


            if (event_log)
                event_log->push_back("Updated parameters: " +
                                     to_string(dlib::rowm(parameters, dlib::range(number_of_inequality_constraints,
                                                                                  number_of_inequality_constraints +
                                                                                  number_of_parameters - 1)),
                                               true));
            cost_vector = functor(dlib::rowm(parameters, dlib::range(number_of_inequality_constraints,
                                                                     number_of_inequality_constraints +
                                                                     number_of_parameters - 1))); // ob

            t(0) = (objective_function_value - cost_vector(0)) / (std::max(std::abs(cost_vector(0)), 1.0));
            objective_function_value = cost_vector(0);

            if (number_of_constraints != 0) {
                constraints = dlib::rowm(cost_vector, dlib::range(1, number_of_constraints));

                if (number_of_inequality_constraints != 0) {

                    if (dlib::min(
                            dlib::rowm(constraints,
                                       dlib::range(number_of_equality_constraints, number_of_constraints - 1)) -
                            dlib::subm(parameter_bounds, dlib::rectangle(0, 0, 0, number_of_inequality_constraints - 1))
                    ) > 0.0 &&
                        dlib::min(
                                dlib::subm(parameter_bounds,
                                           dlib::rectangle(1, 0, 1, number_of_inequality_constraints - 1)) -
                                dlib::rowm(constraints,
                                           dlib::range(number_of_equality_constraints, number_of_constraints - 1))
                        ) > 0.0) {
                        dlib::set_rowm(parameters, dlib::range(0, number_of_inequality_constraints - 1)) =
                                dlib::rowm(constraints,
                                           dlib::range(number_of_equality_constraints, number_of_constraints - 1));
                    }
                    dlib::set_rowm(constraints,
                                   dlib::range(number_of_equality_constraints, number_of_constraints - 1)) =
                            dlib::rowm(constraints,
                                       dlib::range(number_of_equality_constraints, number_of_constraints - 1)) -
                            dlib::rowm(parameters, dlib::range(0, number_of_inequality_constraints - 1));
                }
                t(2) = euclidean_norm(constraints);

                if (t(2) < 10.0 * tolerance) {
                    rho = 0.0;
                    mu = std::min(mu, tolerance);
                }

                if (t(2) < 5.0 * t(1)) {
                    rho /= 5.0;
                } else if (t(2) > 10.0 * t(1)) {
                    rho = 5.0 * std::max(rho, std::sqrt(tolerance));
                }

                if (std::max(tolerance + t(0), t(1) - t(2)) <= 0.0) {
                    lagrangian_multipliers = 0;
                    hessian_matrix = dlib::diagm(dlib::diag(hessian_matrix));
                }
                t(1) = t(2);
            }


            if (std::hypot(t(0), t(1)) <= tolerance) {
                // Tolerance reached, stop procedure
                maximum_major_iterations = iteration;
            }
            objective_function_value_history.push_back(objective_function_value);

        }

        // Save result
        dlib::matrix<double, 0, 1> optimal_parameters =
                dlib::rowm(parameters, dlib::range(
                        number_of_inequality_constraints,
                        number_of_inequality_constraints + number_of_parameters - 1)
                );

        bool converged = false;
        if (std::hypot(t(0), t(1)) <= tolerance) {
            if (event_log) {
                // Reached tolerance
                event_log->push_back("Reached requested tolerance in " + std::to_string(iteration) + " iterations.");
            }
            converged = true;
        } else {
            if (event_log) {
                event_log->push_back("Exiting after maximum number of iterations. Tolerance not reached.");
            }
        }

        if (event_log) {
            event_log->push_back("Final result:" + cppsolnp::to_string(optimal_parameters, false));
        }

        SolveResult result(objective_function_value, optimal_parameters, converged, hessian_matrix);

        return result;
    }


    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    SolveResult solnp(
            functor_model functor,
            parameter_input &parameter_data,
            const inequality_constraint_vectors &inequality_constraint_data,
            // Optional input
            const std::shared_ptr<std::vector<std::string>> &event_log = nullptr,
            double rho = 1.0, //penalty parameter
            int maximum_major_iterations = 400,
            int maximum_minor_iterations = 800,
            const double &delta = 1e-7, // Step size in forward difference evaluation
            const double &tolerance = 1e-8
    ) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<parameter_input>::value);
        COMPILE_TIME_ASSERT(dlib::is_matrix<inequality_constraint_vectors>::value);
        COMPILE_TIME_ASSERT(parameter_input::NC <= 3 && parameter_input::NC >= 0);
        COMPILE_TIME_ASSERT(inequality_constraint_vectors::NC <= 3);

        dlib::matrix<double> hessian_matrix = dlib::identity_matrix<double>(
                parameter_data.nr() + inequality_constraint_data.nr());
        return solnp(functor, parameter_data, inequality_constraint_data, hessian_matrix, event_log, rho,
                     maximum_major_iterations, maximum_minor_iterations, delta, tolerance);
    }

    template<
            typename functor_model,
            typename parameter_input>
    SolveResult solnp(
            functor_model functor,
            parameter_input &parameter_data,
            // Optional input
            const std::shared_ptr<std::vector<std::string>> &event_log = nullptr,
            double rho = 1.0, //penalty parameter
            int maximum_major_iterations = 400,
            int maximum_minor_iterations = 800,
            const double &delta = 1e-7, // Step size in forward difference evaluation
            const double &tolerance = 1e-8
    ) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<parameter_input>::value);
        COMPILE_TIME_ASSERT(parameter_input::NC <= 3 && parameter_input::NC >= 0);

        dlib::matrix<double, 0, 0> inequality_bounds;
        return solnp(functor, parameter_data, inequality_bounds, event_log, rho, maximum_major_iterations,
                     maximum_minor_iterations, delta, tolerance);
    }

}; //namespace cppsolnp

#endif // CPP_SOLNP_SOLNP_HPP
