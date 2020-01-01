#ifndef CPP_SOLNP_SOLVER_HPP
#define CPP_SOLNP_SOLVER_HPP

#include "solnp.hpp"

/*
 * This file holds the entry functions for using CPP SOLNP
 */

namespace cppsolnp {


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