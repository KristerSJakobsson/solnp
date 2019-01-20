//
// Created by krister on 2019/01/20.
//

#ifndef CPP_SOLNP_SUBNP_HPP
#define CPP_SOLNP_SUBNP_HPP

#include "stdafx.h"

#include "utils.hpp"

namespace cppsolnp {


    template<typename functor_model>
    struct subnp {
    public:
        subnp(functor_model &objective_function,
              int number_of_parameter_data, // op(7)
              int number_of_equality_constraints, // op(5)
              int number_of_inequality_constraints,
              const std::pair<bool, bool> &lagrangian_parameters_bounded,//op(6))
              const log_list_ptr &event_log) :
                objective_function_(objective_function),
                number_of_parameters_(number_of_parameter_data), // op(7)
                number_of_equality_constraints_(number_of_equality_constraints),
                number_of_inequality_constraints_(number_of_inequality_constraints),
                number_of_total_constraints_(number_of_inequality_constraints + number_of_equality_constraints),
                number_of_parameters_and_inequality_constraints_(
                        number_of_parameter_data + number_of_inequality_constraints),
                event_log_(event_log), lagrangian_parameters_bounded_(lagrangian_parameters_bounded) {
            /* Here we put subnp initialization data, like function declarations or constant parameter_data. */
        }

        void operator()(dlib::matrix<double> &parameter, // p0
                        dlib::matrix<double> parameter_bounds,
                        dlib::matrix<double, 0, 1> &lagrangian_multipliers,
                        dlib::matrix<double> cost_vector,// ob( )
                        dlib::matrix<double> &hessian, //hessian
                        double &mu,
                        double rho, // op(1)
                        int max_iterations, // op(2)
                        const double delta = 1e-5, // op(3)
                        const double tolerance = 1e-7 // op(4)
        ) {


            /* For debugging*/
            std::string debug_parameter;
            std::string debug_parameter0;
            std::string debug_parameter_bounds;
            std::string debug_lagrangian_multipliers;
            std::string debug_lagrangian_multipliers0;
            std::string debug_hessian;
            std::string debug_cost_vector;
            std::string debug_cost_vector1;
            std::string debug_cost_vector2;
            std::string debug_cost_vector3;
            std::string debug_scale;
            std::string debug_a;
            std::string debug_b;
            std::string debug_c;
            std::string debug_constraints;
            std::string debug_gradient;
            std::string debug_dx;
            std::string debug_gap;
            std::string debug_temporary_vector;
            std::string debug_temporary_gradient;
            std::string debug_modified_cost_vector;
            std::string debug_temporary_parameter;
            std::string debug_u;
            std::string debug_cholesky;
            std::string debug_pt;

            /* Here we put the subnp contents.*/
            alpha_ = dlib::zeros_matrix<double>(3, 1);
            positive_change_ = true;

            dlib::matrix<double> parameter0 = parameter;
            dlib::matrix<double> lagrangian_multipliers0 = lagrangian_multipliers;

            // Debug
            debug_parameter = to_string(parameter);
            debug_parameter0 = to_string(parameter0);
            debug_lagrangian_multipliers = to_string(lagrangian_multipliers);
            debug_lagrangian_multipliers0 = to_string(lagrangian_multipliers0);

            /* Calculate scale for cost, equality constraints,
            inequality constraints and parameter_data.*/
            dlib::matrix<double> scale(1 + number_of_equality_constraints_, 1);
            if (number_of_equality_constraints_ != 0) {
                scale = dlib::join_cols(dlib::mat(cost_vector(0)),
                                        dlib::ones_matrix<double>(number_of_equality_constraints_, 1) *
                                        infinity_norm(dlib::rowm(cost_vector,
                                                                 dlib::range(1, number_of_equality_constraints_)))
                );


                //scale(0) = cost_vector(0);
                //dlib::set_rowm(scale, dlib::range(1, number_of_equality_constraints_)) =
                //	dlib::ones_matrix<double>(number_of_equality_constraints_, 1) *
                //	infinity_norm(dlib::rowm(cost_vector, dlib::range(1, number_of_equality_constraints_)));
            } else {
                scale = dlib::mat(1.0);
            }


            if (lagrangian_parameters_bounded_.second == false) {
                scale = dlib::join_cols(scale, parameter0);
            } else {
                scale = dlib::join_cols(scale, dlib::ones_matrix<double>(parameter0.nr(), parameter0.nc()));
            }

            scale = elementwise_min(elementwise_max(dlib::abs(scale), tolerance), 1.0 / tolerance);

            // Debug
            debug_scale = to_string(scale);


            /*
            scale the cost, the equality constraints, the inequality constraints,
            the parameters (inequality parameters AND actual parameters),
            and the parameter bounds if there are any
            Also make sure the parameters are no larger than (1-tol) times their bounds
            */

            cost_vector = pointwise_divide(cost_vector,
                                           dlib::rowm(scale, dlib::range(0, number_of_total_constraints_)));
            parameter0 = pointwise_divide(parameter0,
                                          dlib::rowm(scale, dlib::range(number_of_equality_constraints_ + 1,
                                                                        number_of_total_constraints_ +
                                                                        number_of_parameters_)));

            // Debug
            debug_cost_vector = to_string(cost_vector);
            debug_parameter0 = to_string(parameter0);

            int mm;
            if (lagrangian_parameters_bounded_.second == true) {
                if (lagrangian_parameters_bounded_.first == false) {
                    mm = number_of_inequality_constraints_;
                } else {
                    mm = number_of_parameters_and_inequality_constraints_;
                }
                /* Scale parameter bounds */
                parameter_bounds = pointwise_divide(parameter_bounds,
                                                    dlib::join_rows(
                                                            dlib::rowm(scale,
                                                                       dlib::range(number_of_equality_constraints_ + 1,
                                                                                   number_of_equality_constraints_ +
                                                                                   mm)),
                                                            dlib::rowm(scale,
                                                                       dlib::range(number_of_equality_constraints_ + 1,
                                                                                   number_of_equality_constraints_ +
                                                                                   mm))
                                                    )
                );
            }

            // Debug
            debug_parameter_bounds = to_string(parameter_bounds);

            if (number_of_total_constraints_ != 0) {
                lagrangian_multipliers0 = (dlib::pointwise_multiply(
                        dlib::rowm(scale, dlib::range(1, number_of_total_constraints_)),
                        lagrangian_multipliers0)) / scale(0);

            }


            hessian = dlib::pointwise_multiply(
                    hessian,
                    dlib::rowm(scale, dlib::range(number_of_equality_constraints_ + 1,
                                                  number_of_total_constraints_ + number_of_parameters_)) *
                    dlib::trans(dlib::rowm(scale, dlib::range(number_of_equality_constraints_ + 1,
                                                              number_of_total_constraints_ + number_of_parameters_)))) /
                      scale(0);

            // Debug
            debug_lagrangian_multipliers0 = to_string(lagrangian_multipliers0);
            debug_hessian = to_string(hessian);

            double object_function_value = cost_vector(0);
            dlib::matrix<double> a(number_of_equality_constraints_ + number_of_inequality_constraints_,
                                   number_of_inequality_constraints_ + number_of_parameters_);

            //a = dlib::join_cols(
            //	dlib::zeros_matrix<double>(number_of_equality_constraints_, number_of_inequality_constraints_),
            //	-1 * dlib::identity_matrix<double>(number_of_inequality_constraints_)
            //);

            dlib::set_colm(a, dlib::range(0, number_of_inequality_constraints_ - 1)) =
                    dlib::join_cols(
                            dlib::zeros_matrix<double>(
                                    number_of_equality_constraints_,
                                    number_of_inequality_constraints_
                            ),
                            -1 * dlib::identity_matrix<double>(number_of_inequality_constraints_)
                    );

            // Debug
            debug_a = to_string(a);


            /*
            Automatic Differentiation
            */
            dlib::matrix<double> gradient;
            gradient = dlib::zeros_matrix<double>(number_of_parameters_and_inequality_constraints_, 1);
            dlib::matrix<double> b;

            dlib::matrix<double> constraints;


            if (number_of_total_constraints_ != 0) {
                constraints = dlib::rowm(cost_vector, dlib::range(1, number_of_total_constraints_));

                // Debug
                debug_constraints = to_string(constraints);

                for (int i = 0; i < number_of_parameters_; ++i) {
                    parameter0(number_of_inequality_constraints_ + i) =
                            parameter0(number_of_inequality_constraints_ + i) + delta;
                    cost_vector =
                            pointwise_divide(
                                    objective_function_(
                                            dlib::pointwise_multiply(
                                                    dlib::rowm(parameter0,
                                                               dlib::range(number_of_inequality_constraints_,
                                                                           number_of_parameters_and_inequality_constraints_ -
                                                                           1)),
                                                    dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                                                  number_of_total_constraints_ +
                                                                                  number_of_parameters_))
                                            )
                                    ),
                                    dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                            );
                    gradient(number_of_inequality_constraints_ + i) =
                            (cost_vector(0) - object_function_value) / delta;

                    dlib::set_colm(a, number_of_inequality_constraints_ + i) =
                            (dlib::rowm(cost_vector, dlib::range(1, number_of_total_constraints_)) - constraints) /
                            delta;


                    parameter0(number_of_inequality_constraints_ + i) =
                            parameter0(number_of_inequality_constraints_ + i) - delta;

                    // Debug
                    debug_parameter0 = to_string(parameter0);
                    debug_cost_vector = to_string(cost_vector);
                    debug_gradient = to_string(gradient);
                    debug_a = to_string(a);

                }
                if (number_of_inequality_constraints_ != 0) {
                    dlib::set_rowm(constraints, dlib::range(number_of_equality_constraints_,
                                                            number_of_equality_constraints_ +
                                                            number_of_inequality_constraints_ - 1)) =
                            dlib::rowm(constraints, dlib::range(number_of_equality_constraints_,
                                                                number_of_equality_constraints_ +
                                                                number_of_inequality_constraints_ - 1)) -
                            dlib::rowm(parameter0, dlib::range(0, number_of_inequality_constraints_ - 1));

                    // Debug
                    debug_constraints = to_string(constraints);
                }
                if (conditional_number(a) > 1 / std::numeric_limits<double>::epsilon()) {
                    event_log_->push_back(
                            "Warning: Redundant constraints were detected. Poor intermediate results may result.");


                }
                // Debug
                double cond_numb = conditional_number(a);
                b = a * parameter0 - constraints;

                // Debug
                debug_b = to_string(b);
            }


            dlib::matrix<double> c(1, number_of_parameters_and_inequality_constraints_ + 1);
            dlib::matrix<double> dx(number_of_parameters_and_inequality_constraints_ + 1, 1);
            double go;
            int minor_iteration;
            int major_iteration;
            dlib::matrix<double> gap(parameter_bounds.nr(), 2);
            if (number_of_total_constraints_ != 0) {
                positive_change_ = false; // ch = -1
                alpha_(0) = tolerance - dlib::max(dlib::abs(constraints));

                // Debug alpha

                if (alpha_(0) <= 0) {
                    positive_change_ = true;
                    if (lagrangian_parameters_bounded_.second == false) {
                        /*pseudo inverse?*/
                        dlib::qr_decomposition <dlib::matrix<double>> qr_temp(a * trans(a));

                        parameter0 = parameter0 - trans(a) * (qr_temp.solve(constraints));
                        alpha_(0) = 1;

                        // Debug
                        debug_parameter0 = to_string(parameter0);
                    }
                }

                if (alpha_(0) <= 0) {
                    parameter0 = dlib::join_cols(parameter0, dlib::mat(1.0));

                    a = dlib::join_rows(a, -1 * constraints);
                    c = dlib::zeros_matrix<double>(1, number_of_parameters_and_inequality_constraints_ + 1);
                    c(number_of_parameters_and_inequality_constraints_) = 1.0;
                    dx = dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_ + 1, 1);
                    go = 1.0;
                    minor_iteration = 0;

                    // Debug
                    debug_a = to_string(a);
                    debug_c = to_string(c);
                    debug_dx = to_string(dx);

                    while (go >= tolerance) {
                        ++minor_iteration;
                        gap = dlib::join_rows(
                                dlib::rowm(parameter0, dlib::range(0, mm - 1)) -
                                dlib::colm(parameter_bounds, 0),
                                dlib::colm(parameter_bounds, 1) -
                                dlib::rowm(parameter0, dlib::range(0, mm - 1))
                        );

                        // Debug
                        debug_gap = to_string(gap);

                        two_vector_sort(gap);

                        // Debug
                        debug_gap = to_string(gap);

                        dlib::set_rowm(dx, dlib::range(0, mm - 1)) = dlib::colm(gap, 0);
                        dx(number_of_parameters_and_inequality_constraints_) =
                                parameter0(number_of_parameters_and_inequality_constraints_);

                        // Debug
                        debug_dx = to_string(dx);

                        if (lagrangian_parameters_bounded_.first == false) {
                            dlib::set_rowm(dx, dlib::range(mm, number_of_parameters_and_inequality_constraints_ - 1)) =
                                    std::max(dlib::max(dlib::rowm(dx, dlib::range(0, mm - 1))), 100.0) *
                                    dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_ - mm, 1);
                            // Debug
                            debug_dx = to_string(dx);
                        }

                        // TODO: Ought to be better way than to use the Pseudo Inverse.


#ifdef QR
                        dlib::qr_decomposition<dlib::matrix<double>> qr(dlib::trans(a * dlib::diagm(dx)));
                            lagrangian_multipliers = qr.solve(dlib::pointwise_multiply(dx, dlib::trans(c)));

#endif

#ifdef SVD
                        lagrangian_multipliers = dlib::pinv(dlib::trans(a*dlib::diagm(dx)))*dlib::pointwise_multiply(dx, dlib::trans(c));
#endif
                        // Debug
                        debug_lagrangian_multipliers = to_string(lagrangian_multipliers);

                        dlib::matrix<double, 0, 1> temporary_vector(dx.nr()); // v

                        temporary_vector = dlib::pointwise_multiply(dx,
                                                                    dlib::pointwise_multiply(dx,
                                                                                             dlib::trans(c) -
                                                                                             dlib::trans(a) *
                                                                                             lagrangian_multipliers));
                        // Debug
                        debug_temporary_vector = to_string(temporary_vector);

                        if (temporary_vector(number_of_parameters_and_inequality_constraints_) > 0) {
                            double temporary_scalar; // z
                            temporary_scalar = parameter0(number_of_parameters_and_inequality_constraints_) /
                                               temporary_vector(number_of_parameters_and_inequality_constraints_);

                            // Debug temporary scalar
                            for (int k = 0; k < mm; k++) {
                                if (temporary_vector(k) < 0) {
                                    temporary_scalar = std::min(temporary_scalar,
                                                                -1 * (parameter_bounds(k, 1) - parameter0(k)) /
                                                                temporary_vector(k));
                                } else if (temporary_vector(k) > 0) {
                                    temporary_scalar = std::min(temporary_scalar,
                                                                (parameter0(k) - parameter_bounds(k, 0)) /
                                                                temporary_vector(k));
                                }
                            }

                            // Debug temporary scalar

                            if (temporary_scalar >= parameter0(number_of_parameters_and_inequality_constraints_) /
                                                    temporary_vector(
                                                            number_of_parameters_and_inequality_constraints_)) {
                                parameter0 = parameter0 - temporary_scalar * temporary_vector;
                            } else {
                                parameter0 = parameter0 - 0.9 * temporary_scalar * temporary_vector;
                            }
                            // Debug
                            debug_parameter0 = to_string(parameter0);


                            go = parameter0(number_of_parameters_and_inequality_constraints_);
                            if (minor_iteration >= 10) {
                                go = 0.0;
                            }
                        } else {
                            go = 0.0;
                            minor_iteration = 10;
                        }
                    }
                    if (minor_iteration >= 10) {
                        event_log_->push_back(
                                "Warning: The linearized prblem has no feasible solution. The problem may not be feasible.");
                    }
                    a = dlib::colm(a, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));
                    b = a *
                        dlib::rowm(parameter0, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));

                    // Debug
                    debug_a = to_string(a);
                    debug_b = to_string(b);


                }
            }

            parameter = dlib::rowm(parameter0, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));
            lagrangian_multipliers = 0;
            // Debug
            debug_parameter = to_string(parameter);
            debug_lagrangian_multipliers = to_string(lagrangian_multipliers);

            if (positive_change_ == true) {
                cost_vector = pointwise_divide(
                        objective_function_(
                                dlib::pointwise_multiply(
                                        dlib::rowm(parameter, dlib::range(number_of_inequality_constraints_,
                                                                          number_of_parameters_and_inequality_constraints_ -
                                                                          1)),
                                        dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                                      number_of_total_constraints_ +
                                                                      number_of_parameters_)))),
                        dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                );

                // Debug
                debug_cost_vector = to_string(cost_vector);
            }

            object_function_value = cost_vector(0);

            if (number_of_inequality_constraints_ != 0) {
                dlib::set_rowm(cost_vector,
                               dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                        dlib::rowm(cost_vector,
                                   dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) -
                        dlib::rowm(parameter, dlib::range(0, number_of_inequality_constraints_ - 1));
                // Debug
                debug_cost_vector = to_string(cost_vector);
            }

            if (number_of_total_constraints_ != 0) {
                dlib::set_rowm(cost_vector, dlib::range(1, number_of_total_constraints_)) =
                        dlib::rowm(cost_vector, dlib::range(1, number_of_total_constraints_)) -
                        a * parameter + b;


                object_function_value = cost_vector(0) - dlib::trans(lagrangian_multipliers0) * dlib::rowm(cost_vector,
                                                                                                           dlib::range(
                                                                                                                   1,
                                                                                                                   number_of_total_constraints_)) +
                                        rho * std::pow<double>(euclidean_norm(
                                                dlib::rowm(cost_vector, dlib::range(1, number_of_total_constraints_))),
                                                               2);

                // Debug
                debug_cost_vector = to_string(cost_vector);

            }
            minor_iteration = 0;

            dlib::matrix<double> modified_cost_vector;
            dlib::matrix<double> temporary_gradient, temporary_parameter;
            temporary_gradient = dlib::zeros_matrix<double>(parameter.nr(), parameter.nc());
            temporary_parameter = dlib::zeros_matrix<double>(parameter.nr(), parameter.nc());
            double modified_object_function_value;
            double reduction;
            while (minor_iteration < max_iterations) {
                ++minor_iteration;
                if (positive_change_ == true) {

                    // H�r �r n�got fel! Loopen �ndrar g konstigt
                    for (auto i = 0; i < number_of_parameters_; ++i) {
                        parameter(number_of_inequality_constraints_ + i) =
                                parameter(number_of_inequality_constraints_ + i) + delta;

                        modified_cost_vector = pointwise_divide(
                                objective_function_(
                                        dlib::pointwise_multiply(
                                                dlib::rowm(parameter, dlib::range(number_of_inequality_constraints_,
                                                                                  number_of_parameters_and_inequality_constraints_ -
                                                                                  1)),
                                                dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                                              number_of_total_constraints_ +
                                                                              number_of_parameters_)))),
                                dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                        );

                        // Debug
                        debug_parameter = to_string(parameter);
                        debug_modified_cost_vector = to_string(modified_cost_vector);

                        if (number_of_inequality_constraints_ != 0) {
                            dlib::set_rowm(modified_cost_vector, dlib::range(number_of_equality_constraints_ + 1,
                                                                             number_of_total_constraints_)) =
                                    dlib::rowm(modified_cost_vector, dlib::range(number_of_equality_constraints_ + 1,
                                                                                 number_of_total_constraints_)) -
                                    dlib::rowm(parameter, dlib::range(0, number_of_inequality_constraints_ - 1));

                            // Debug
                            debug_modified_cost_vector = to_string(modified_cost_vector);
                        }

                        if (number_of_total_constraints_ != 0) {
                            dlib::set_rowm(modified_cost_vector, dlib::range(1, number_of_total_constraints_)) =
                                    dlib::rowm(modified_cost_vector, dlib::range(1, number_of_total_constraints_)) -
                                    a * parameter + b;

                            // Debug
                            debug_modified_cost_vector = to_string(modified_cost_vector);
                            double test = euclidean_norm(
                                    dlib::rowm(modified_cost_vector, dlib::range(1, number_of_total_constraints_)));

                            modified_object_function_value = modified_cost_vector(0) -
                                                             dlib::trans(lagrangian_multipliers0) *
                                                             dlib::rowm(modified_cost_vector,
                                                                        dlib::range(1, number_of_total_constraints_)) +
                                                             rho * std::pow<double>(euclidean_norm(
                                                                     dlib::rowm(modified_cost_vector, dlib::range(1,
                                                                                                                  number_of_total_constraints_))),
                                                                                    2);


                        }

                        gradient(number_of_inequality_constraints_ + i) =
                                (modified_object_function_value - object_function_value) / delta;
                        parameter(number_of_inequality_constraints_ + i) =
                                parameter(number_of_inequality_constraints_ + i) - delta;

                        // Debug
                        debug_gradient = to_string(gradient);
                        debug_parameter = to_string(parameter);
                    }
                    if (number_of_inequality_constraints_ != 0) {
                        dlib::set_rowm(gradient, dlib::range(0, number_of_inequality_constraints_ - 1)) =
                                dlib::zeros_matrix<double>(number_of_inequality_constraints_, 1);

                        // Debug
                        debug_gradient = to_string(gradient);
                    }
                }

                if (minor_iteration > 1) {
                    temporary_gradient = gradient - temporary_gradient;
                    temporary_parameter = parameter - temporary_parameter;

                    // Debug
                    debug_temporary_gradient = to_string(gradient);
                    debug_temporary_parameter = to_string(gradient);

                    double sc[2];
                    sc[0] = dlib::trans(temporary_parameter) * hessian * temporary_parameter;
                    sc[1] = dlib::trans(temporary_parameter) * temporary_gradient;
                    if (sc[0] * sc[1] > 0) {
                        temporary_parameter = hessian * temporary_parameter;
                        hessian = hessian - (temporary_parameter * dlib::trans(temporary_parameter)) / sc[0] +
                                  (temporary_gradient * dlib::trans(temporary_gradient)) / sc[1];

                    }
                    // Debug
                    debug_temporary_parameter = to_string(temporary_parameter);
                    debug_hessian = to_string(hessian);
                }

                dx = 0.01 * dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_, 1);
                // Debug
                debug_dx = to_string(dx);

                if (lagrangian_parameters_bounded_.second == true) {
                    gap = dlib::join_rows(
                            dlib::rowm(parameter, dlib::range(0, mm - 1)) -
                            dlib::colm(parameter_bounds, 0),
                            dlib::colm(parameter_bounds, 1) -
                            dlib::rowm(parameter, dlib::range(0, mm - 1))
                    );

                    // Debug
                    debug_gap = to_string(gap);

                    two_vector_sort(gap);

                    // Debug
                    debug_gap = to_string(gap);

                    gap = dlib::colm(gap, 0) +
                          std::sqrt(std::numeric_limits<double>::epsilon()) * dlib::ones_matrix<double>(mm, 1);

                    // Debug
                    debug_gap = to_string(gap);

                    dlib::set_rowm(dx, dlib::range(0, mm - 1)) = pointwise_divide(dlib::ones_matrix<double>(mm, 1),
                                                                                  gap);

                    // Debug
                    debug_dx = to_string(dx);

                    if (lagrangian_parameters_bounded_.first == false) {
                        dlib::set_rowm(dx, dlib::range(mm, number_of_parameters_and_inequality_constraints_ - 1)) =
                                dlib::min(dlib::join_cols(dlib::rowm(dx, dlib::range(0, mm - 1)), dlib::mat(0.01))) *
                                dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_ - mm, 1);

                        // Debug
                        debug_dx = to_string(dx);

                    }
                }
                go = -1.0;
                mu = mu / 10.0;
                dlib::matrix<double> u;
                dlib::matrix<double> cholesky;
                while (go <= 0) {

                    //if (chole.is_spd() == false)
                    //{
                    //	//throw std::runtime_error("Fail");
                    //	event_log_->push_back("Warning: Cholesky decomposition failed. Return best found value.");
                    //	//cholesky = cholesky_last;
                    //	break;
                    //}

                    cholesky = dlib::trans(dlib::chol(
                            dlib::make_symmetric(hessian) + mu * dlib::diagm(dlib::pointwise_multiply(dx, dx))));
                    //dlib::trans(dlib::chol(hessian + mu*dlib::diagm(dlib::pointwise_multiply(dx, dx))));

                    // Debug
                    debug_cholesky = to_string(cholesky);

#ifdef CHOL_UPP
                    cholesky = dlib::inv_upper_triangular(cholesky);
#endif

#ifdef CHOL_WHOLE
                    cholesky = dlib::inv(cholesky);
#endif
#ifdef CHOL_QR
                    dlib::qr_decomposition<dlib::matrix<double>> chol_qr;

                        cholesky = dlib::inv(cholesky);
#endif

                    // Debug
                    debug_cholesky = to_string(cholesky);

                    temporary_gradient = dlib::trans(cholesky) * gradient;

                    // Debug
                    debug_temporary_gradient = to_string(temporary_gradient);

                    if (number_of_total_constraints_ == 0) {
                        u = -1 * cholesky * temporary_gradient;

                        // Debug
                        debug_u = to_string(u);

                    } else {
                        // We solve the equation system using QR factorization
#ifdef QR
                        dlib::qr_decomposition<dlib::matrix<double>> qr(trans(cholesky) * dlib::trans(a));
                            lagrangian_multipliers = qr.solve(temporary_gradient);
#endif

#ifdef SVD
                        lagrangian_multipliers = dlib::pinv(trans(cholesky)*dlib::trans(a))*temporary_gradient;
#endif


                        u = -1 * cholesky *
                            (temporary_gradient - (dlib::trans(cholesky) * dlib::trans(a)) * lagrangian_multipliers);
                        // Debug
                        debug_lagrangian_multipliers = to_string(lagrangian_multipliers);
                        debug_u = to_string(u);
                    }

                    parameter0 = dlib::rowm(u, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1)) +
                                 parameter;

                    // Debug
                    debug_parameter0 = to_string(parameter0);

                    if (lagrangian_parameters_bounded_.second == false) {
                        go = 1.0;
                    } else {
                        go = dlib::min(
                                dlib::join_cols(
                                        dlib::rowm(parameter0, dlib::range(0, mm - 1)) -
                                        dlib::colm(parameter_bounds, 0),
                                        dlib::colm(parameter_bounds, 1) - dlib::rowm(parameter0, dlib::range(0, mm - 1))
                                )
                        );
                        mu *= 3.0;
                    }
                    // Debug go

                }
                dlib::matrix<double> cost_vector1, cost_vector2, cost_vector3;
                dlib::matrix<double> pt(parameter.nr(), 3);
                dlib::matrix<double, 3, 1> sob;

                alpha_(0) = 0;

                cost_vector1 = cost_vector;
                cost_vector2 = cost_vector;
                sob(0) = object_function_value;
                sob(1) = object_function_value;
                dlib::set_colm(pt, dlib::range(0, 1)) = dlib::join_rows(parameter, parameter);
                alpha_(2) = 1.0;
                dlib::set_colm(pt, 2) = parameter0;
                cost_vector3 = pointwise_divide(
                        objective_function_(
                                dlib::pointwise_multiply(
                                        dlib::subm(pt, dlib::rectangle(2, number_of_inequality_constraints_, 2,
                                                                       number_of_parameters_and_inequality_constraints_ -
                                                                       1)),
                                        dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                                      number_of_total_constraints_ +
                                                                      number_of_parameters_)))),
                        dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                );
                sob(2) = cost_vector3(0);

                // Debug
                debug_cost_vector1 = to_string(cost_vector1);
                debug_cost_vector2 = to_string(cost_vector2);
                debug_cost_vector3 = to_string(cost_vector3);
                debug_pt = to_string(pt);


                if (number_of_inequality_constraints_ != 0) {
                    dlib::set_rowm(cost_vector3,
                                   dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                            dlib::rowm(cost_vector3,
                                       dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) -
                            dlib::subm(pt, dlib::rectangle(2, 0, 2, number_of_inequality_constraints_ - 1));

                    // Debug
                    debug_cost_vector3 = to_string(cost_vector3);
                }
                if (number_of_total_constraints_ != 0) {
                    dlib::set_rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) =
                            dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) -
                            a * dlib::colm(pt, 2) + b;
                    sob(2) = cost_vector3(0) - dlib::trans(lagrangian_multipliers0) *
                                               dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) +
                             rho * std::pow<double>(euclidean_norm(
                                     dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_))), 2);
                    // Debug
                    debug_cost_vector3 = to_string(cost_vector3);
                }

                go = 1.0;
                double obm, obn;
                while (go > tolerance) {
                    alpha_(1) = 0.5 * (alpha_(0) + alpha_(2));
                    dlib::set_colm(pt, 1) = (1 - alpha_(1)) * parameter + alpha_(1) * parameter0;


                    cost_vector2 = pointwise_divide(
                            objective_function_(
                                    dlib::pointwise_multiply(
                                            dlib::subm(pt, dlib::rectangle(1, number_of_inequality_constraints_, 1,
                                                                           number_of_parameters_and_inequality_constraints_ -
                                                                           1)),
                                            dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                                          number_of_total_constraints_ +
                                                                          number_of_parameters_)))),
                            dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                    );
                    sob(1) = cost_vector2(0);

                    // Debug
                    debug_pt = to_string(pt);
                    debug_cost_vector2 = to_string(cost_vector2);

                    if (number_of_inequality_constraints_ != 0) {
                        dlib::set_rowm(cost_vector2,
                                       dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                                dlib::rowm(cost_vector2, dlib::range(number_of_equality_constraints_ + 1,
                                                                     number_of_total_constraints_)) -
                                dlib::subm(pt, dlib::rectangle(1, 0, 1, number_of_inequality_constraints_ - 1));
                        // Debug
                        debug_cost_vector2 = to_string(cost_vector2);

                    }
                    if (number_of_total_constraints_ != 0) {
                        dlib::set_rowm(cost_vector2, dlib::range(1, number_of_total_constraints_)) =
                                dlib::rowm(cost_vector2, dlib::range(1, number_of_total_constraints_)) -
                                a * dlib::colm(pt, 1) + b;
                        sob(1) = cost_vector2(0) - dlib::trans(lagrangian_multipliers0) *
                                                   dlib::rowm(cost_vector2,
                                                              dlib::range(1, number_of_total_constraints_)) +
                                 rho * std::pow<double>(euclidean_norm(
                                         dlib::rowm(cost_vector2, dlib::range(1, number_of_total_constraints_))), 2);
                        // Debug
                        debug_cost_vector2 = to_string(cost_vector2);
                    }
                    obm = dlib::max(sob);
                    if (obm < object_function_value) {
                        obn = dlib::min(sob);
                        go = tolerance * (obm - obn) / (object_function_value - obm);

                        // Debug go & obn

                    }

                    if (sob(1) >= sob(0)) {
                        sob(2) = sob(1);
                        cost_vector3 = cost_vector2;
                        alpha_(2) = alpha_(1);
                        dlib::set_colm(pt, 2) = dlib::colm(pt, 1);

                        // Debug
                        debug_pt = to_string(pt);
                        debug_cost_vector3 = to_string(cost_vector3);

                    } else if (sob(0) <= sob(2)) {
                        sob(2) = sob(1);
                        cost_vector3 = cost_vector2;
                        alpha_(2) = alpha_(1);
                        dlib::set_colm(pt, 2) = dlib::colm(pt, 1);

                        // Debug
                        debug_pt = to_string(pt);
                        debug_cost_vector2 = to_string(cost_vector2);
                    } else {
                        sob(0) = sob(1);
                        cost_vector1 = cost_vector2;
                        alpha_(0) = alpha_(1);
                        dlib::set_colm(pt, 0) = dlib::colm(pt, 1);

                        // Debug
                        debug_pt = to_string(pt);
                        debug_cost_vector1 = to_string(cost_vector1);
                    }

                    if (go >= tolerance) {
                        go = alpha_(2) - alpha_(0);
                    }

                    // Debug

                }
                temporary_parameter = parameter;
                temporary_gradient = gradient;
                positive_change_ = true;
                obn = dlib::min(sob);

                // Debug
                debug_temporary_parameter = to_string(temporary_parameter);
                debug_temporary_gradient = to_string(temporary_gradient);

                if (object_function_value < obn) {
                    max_iterations = minor_iteration;

                    // Debug
                }
                reduction = (object_function_value - obn) / (1 + std::abs(object_function_value));
                //Reduction too low? Then we end the loop.
                if (reduction < tolerance) {
                    max_iterations = minor_iteration;
                    // Debug

                }
                if (sob(0) < sob(1)) {
                    object_function_value = sob(0);
                    parameter = dlib::colm(pt, 0);
                    cost_vector = cost_vector1;

                    // Debug
                    debug_parameter = to_string(parameter);
                    debug_cost_vector = to_string(cost_vector);

                } else if (sob(2) < sob(1)) {
                    object_function_value = sob(2);
                    parameter = dlib::colm(pt, 2);
                    cost_vector = cost_vector3;
                    // Debug
                    debug_parameter = to_string(parameter);
                    debug_cost_vector = to_string(cost_vector);
                } else {
                    object_function_value = sob(1);
                    parameter = dlib::colm(pt, 1);
                    cost_vector = cost_vector2;
                    // Debug
                    debug_parameter = to_string(parameter);
                    debug_cost_vector = to_string(cost_vector);
                }
            }
            parameter = dlib::pointwise_multiply(
                    parameter,
                    dlib::rowm(scale, dlib::range(number_of_equality_constraints_ + 1,
                                                  number_of_total_constraints_ + number_of_parameters_))
            );
            // Debug
            debug_parameter = to_string(parameter);


            if (number_of_total_constraints_ != 0) {
                lagrangian_multipliers = scale(0) * pointwise_divide(lagrangian_multipliers,
                                                                     dlib::rowm(scale, dlib::range(1,
                                                                                                   number_of_total_constraints_)));
                // Debug
                debug_lagrangian_multipliers = to_string(lagrangian_multipliers);

            }


            hessian = scale(0) * pointwise_divide(hessian,
                                                  dlib::tmp(dlib::rowm(scale,
                                                                       dlib::range(number_of_equality_constraints_ + 1,
                                                                                   number_of_total_constraints_ +
                                                                                   number_of_parameters_)) *
                                                            dlib::trans(dlib::rowm(scale, dlib::range(
                                                                    number_of_equality_constraints_ + 1,
                                                                    number_of_total_constraints_ +
                                                                    number_of_parameters_))))
            );
            hessian = dlib::make_symmetric(hessian);

            // Debug
            debug_hessian = to_string(hessian);

            if (reduction > tolerance) {
                if (event_log_ != nullptr)
                    event_log_->push_back(
                            "Warning: Minor optimization routine did not converge. You may need to increase the number of minor iterations.");
                // Debug
            }

            std::string debug_return1 = to_string(parameter);
            std::string debug_return2 = to_string(lagrangian_multipliers);
            std::string debug_return3 = to_string(hessian);
            std::string debug_return4 = to_string(scale);
            std::string debug_return5 = std::to_string(mu);

            return;

        }

    private:
        // Constructor variables
        functor_model &objective_function_;
        const int number_of_parameters_;
        const int number_of_equality_constraints_;
        const int number_of_inequality_constraints_;
        const int number_of_total_constraints_;
        const int number_of_parameters_and_inequality_constraints_;
        const std::pair<bool, bool> &lagrangian_parameters_bounded_;
        const log_list_ptr &event_log_;

        // Internal variables
        bool positive_change_; // if true, there have been positive change
        dlib::matrix<double, 3, 1> alpha_;
    };

} // namespace cppsolnp

#endif //CPP_SOLNP_SUBNP_HPP
