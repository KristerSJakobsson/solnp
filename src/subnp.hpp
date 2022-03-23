#ifndef CPP_SOLNP_SUBNP_HPP
#define CPP_SOLNP_SUBNP_HPP

#include "stdafx.h"

#include "utils.hpp"

namespace cppsolnp {


    template<typename functor_model>
    struct subnp {
    public:
        subnp(functor_model &objective_function,
              long number_of_parameter_data,
              long number_of_equality_constraints,
              long number_of_inequality_constraints,
              const std::pair<bool, bool> &lagrangian_parameters_bounded) :
                objective_function_(objective_function),
                number_of_parameters_(number_of_parameter_data),
                number_of_equality_constraints_(number_of_equality_constraints),
                number_of_inequality_constraints_(number_of_inequality_constraints),
                number_of_total_constraints_(number_of_inequality_constraints + number_of_equality_constraints),
                number_of_parameters_and_inequality_constraints_(
                        number_of_parameter_data + number_of_inequality_constraints),
                lagrangian_parameters_bounded_(lagrangian_parameters_bounded) {
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
                        const double tolerance = 1e-7, // op(4)
                        std::shared_ptr<std::vector<std::string>> event_log = nullptr
        ) {


            alpha_ = dlib::zeros_matrix<double>(3, 1);
            bool positive_change = true;

            dlib::matrix<double> parameter0 = parameter;
            dlib::matrix<double> lagrangian_multipliers0 = lagrangian_multipliers;

            /* Calculate scale for cost, equality constraints,
            inequality constraints and parameter_data.*/
            dlib::matrix<double> scale(1 + number_of_equality_constraints_, 1);
            if (number_of_equality_constraints_ != 0) {
                scale = dlib::join_cols(dlib::mat(cost_vector(0)),
                                        dlib::ones_matrix<double>(number_of_equality_constraints_, 1) *
                                        infinity_norm(dlib::rowm(cost_vector,
                                                                 dlib::range(1, number_of_equality_constraints_)))
                );

            } else {
                scale = dlib::mat(1.0);
            }


            if (!lagrangian_parameters_bounded_.second) {
                scale = dlib::join_cols(scale, parameter0);
            } else {
                scale = dlib::join_cols(scale, dlib::ones_matrix<double>(parameter0.nr(), parameter0.nc()));
            }

            scale = elementwise_min(elementwise_max(dlib::abs(scale), tolerance), 1.0 / tolerance);


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

            long mm = 0;
            if (lagrangian_parameters_bounded_.second) {
                if (!lagrangian_parameters_bounded_.first) {
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

            double object_function_value = cost_vector(0);
            dlib::matrix<double> a(number_of_equality_constraints_ + number_of_inequality_constraints_,
                                   number_of_inequality_constraints_ + number_of_parameters_);

            if (number_of_inequality_constraints_ > 0) {
                dlib::set_colm(a, dlib::range(0, number_of_inequality_constraints_ - 1)) =
                        dlib::join_cols(
                                dlib::zeros_matrix<double>(
                                        number_of_equality_constraints_,
                                        number_of_inequality_constraints_
                                ),
                                -1 * dlib::identity_matrix<double>(number_of_inequality_constraints_)
                        );
            }

            /*
            Automatic Differentiation
            */
            dlib::matrix<double> gradient;
            gradient = dlib::zeros_matrix<double>(number_of_parameters_and_inequality_constraints_, 1);
            dlib::matrix<double> b;

            dlib::matrix<double> constraints;


            if (number_of_total_constraints_ != 0) {
                constraints = dlib::rowm(cost_vector, dlib::range(1, number_of_total_constraints_));

                for (auto i = 0; i < number_of_parameters_; ++i) {
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

                }
                if (number_of_inequality_constraints_ != 0) {
                    dlib::set_rowm(constraints, dlib::range(number_of_equality_constraints_,
                                                            number_of_equality_constraints_ +
                                                            number_of_inequality_constraints_ - 1)) =
                            dlib::rowm(constraints, dlib::range(number_of_equality_constraints_,
                                                                number_of_equality_constraints_ +
                                                                number_of_inequality_constraints_ - 1)) -
                            dlib::rowm(parameter0, dlib::range(0, number_of_inequality_constraints_ - 1));

                }

                if (event_log && conditional_number(a) > 1 / std::numeric_limits<double>::epsilon()) {
                    event_log->push_back(
                            "Warning: Redundant constraints were detected. Poor intermediate results may result.");
                }

                b = a * parameter0 - constraints;

            }


            dlib::matrix<double> c(1, number_of_parameters_and_inequality_constraints_ + 1);
            dlib::matrix<double> dx(number_of_parameters_and_inequality_constraints_ + 1, 1);
            double go;
            int minor_iteration = 0;
            dlib::matrix<double> gap(parameter_bounds.nr(), 2);
            if (number_of_total_constraints_ != 0) {
                positive_change = false;
                alpha_(0) = tolerance - dlib::max(dlib::abs(constraints));

                if (alpha_(0) <= 0) {
                    positive_change = true;
                    if (!lagrangian_parameters_bounded_.second) {
                        dlib::qr_decomposition<dlib::matrix<double>> qr_temp(a * trans(a));
                        if (!qr_temp.is_full_rank()) {
                            throw std::domain_error("Encountered Singular matrix when trying to solve. This can happen for example if you have contradicting equality constraints.");
                        }
                        parameter0 = parameter0 - trans(a) * (qr_temp.solve(constraints));
                        alpha_(0) = 1;
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

                    while (go >= tolerance) {
                        ++minor_iteration;
                        gap = dlib::join_rows(
                                dlib::rowm(parameter0, dlib::range(0, mm - 1)) -
                                dlib::colm(parameter_bounds, 0),
                                dlib::colm(parameter_bounds, 1) -
                                dlib::rowm(parameter0, dlib::range(0, mm - 1))
                        );

                        gap = left_vector_min_right_vector_max(gap);

                        dlib::set_rowm(dx, dlib::range(0, mm - 1)) = dlib::colm(gap, 0);
                        dx(number_of_parameters_and_inequality_constraints_) =
                                parameter0(number_of_parameters_and_inequality_constraints_);

                        if (!lagrangian_parameters_bounded_.first) {
                            dlib::set_rowm(dx, dlib::range(mm, number_of_parameters_and_inequality_constraints_ - 1)) =
                                    std::max(dlib::max(dlib::rowm(dx, dlib::range(0, mm - 1))), 100.0) *
                                    dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_ - mm, 1);
                        }


                        dlib::qr_decomposition<dlib::matrix<double>> qr(dlib::trans(a * dlib::diagm(dx)));
                        lagrangian_multipliers = qr.solve(dlib::pointwise_multiply(dx, dlib::trans(c)));

                        dlib::matrix<double, 0, 1> temporary_vector(dx.nr()); // v

                        temporary_vector = dlib::pointwise_multiply(dx,
                                                                    dlib::pointwise_multiply(dx,
                                                                                             dlib::trans(c) -
                                                                                             dlib::trans(a) *
                                                                                             lagrangian_multipliers));

                        if (temporary_vector(number_of_parameters_and_inequality_constraints_) > 0) {
                            double temporary_scalar = parameter0(number_of_parameters_and_inequality_constraints_) /
                                                      temporary_vector(
                                                              number_of_parameters_and_inequality_constraints_);

                            for (auto k = 0L; k < mm; k++) {
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

                            if (temporary_scalar >= parameter0(number_of_parameters_and_inequality_constraints_) /
                                                    temporary_vector(
                                                            number_of_parameters_and_inequality_constraints_)) {
                                parameter0 = parameter0 - temporary_scalar * temporary_vector;
                            } else {
                                parameter0 = parameter0 - 0.9 * temporary_scalar * temporary_vector;
                            }

                            go = parameter0(number_of_parameters_and_inequality_constraints_);
                            if (minor_iteration >= 10) {
                                go = 0.0;
                            }
                        } else {
                            go = 0.0;
                            minor_iteration = 10;
                        }
                    }
                    if (event_log && minor_iteration >= 10) {
                        event_log->push_back(
                                "Warning: The linearized problem has no feasible solution. The problem may not be feasible.");
                    }
                    a = dlib::colm(a, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));
                    b = a *
                        dlib::rowm(parameter0, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));

                }
            }

            parameter = dlib::rowm(parameter0, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1));
            lagrangian_multipliers = 0;

            if (positive_change) {

                dlib::matrix<double, 0, 0> scaled_function_values = dlib::pointwise_multiply(
                        dlib::rowm(parameter, dlib::range(number_of_inequality_constraints_,
                                                          number_of_parameters_and_inequality_constraints_ -
                                                          1)),
                        dlib::rowm(scale, dlib::range(number_of_total_constraints_ + 1,
                                                      number_of_total_constraints_ +
                number_of_parameters_)));

                cost_vector = pointwise_divide(
                        objective_function_(scaled_function_values),
                        dlib::rowm(scale, dlib::range(0, number_of_total_constraints_))
                        );
            }

            object_function_value = cost_vector(0);

            if (number_of_inequality_constraints_ != 0) {
                dlib::set_rowm(cost_vector,
                               dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                        dlib::rowm(cost_vector,
                                   dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) -
                        dlib::rowm(parameter, dlib::range(0, number_of_inequality_constraints_ - 1));
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


            }
            minor_iteration = 0;

            dlib::matrix<double> modified_cost_vector;
            dlib::matrix<double> temporary_gradient, temporary_parameter;
            temporary_gradient = dlib::zeros_matrix<double>(parameter.nr(), parameter.nc());
            temporary_parameter = dlib::zeros_matrix<double>(parameter.nr(), parameter.nc());
            double modified_object_function_value = 0.0;
            double reduction = 0.0;
            while (minor_iteration < max_iterations) {
                ++minor_iteration;
                if (positive_change) {

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

                        if (number_of_inequality_constraints_ != 0) {
                            dlib::set_rowm(modified_cost_vector, dlib::range(number_of_equality_constraints_ + 1,
                                                                             number_of_total_constraints_)) =
                                    dlib::rowm(modified_cost_vector, dlib::range(number_of_equality_constraints_ + 1,
                                                                                 number_of_total_constraints_)) -
                                    dlib::rowm(parameter, dlib::range(0, number_of_inequality_constraints_ - 1));

                        }

                        if (number_of_total_constraints_ != 0) {
                            dlib::set_rowm(modified_cost_vector, dlib::range(1, number_of_total_constraints_)) =
                                    dlib::rowm(modified_cost_vector, dlib::range(1, number_of_total_constraints_)) -
                                    a * parameter + b;

                            modified_object_function_value = modified_cost_vector(0) -
                                                             dlib::trans(lagrangian_multipliers0) *
                                                             dlib::rowm(modified_cost_vector,
                                                                        dlib::range(1, number_of_total_constraints_)) +
                                                             rho * std::pow<double>(euclidean_norm(
                                                                                            dlib::rowm(modified_cost_vector, dlib::range(1,
                                                                                                                                         number_of_total_constraints_))),
                                                                                    2);


                        } else {
                            modified_object_function_value = modified_cost_vector(0);
                        }

                        gradient(number_of_inequality_constraints_ + i) =
                                (modified_object_function_value - object_function_value) / delta;
                        parameter(number_of_inequality_constraints_ + i) =
                                parameter(number_of_inequality_constraints_ + i) - delta;

                    }
                    if (number_of_inequality_constraints_ != 0) {
                        dlib::set_rowm(gradient, dlib::range(0, number_of_inequality_constraints_ - 1)) =
                                dlib::zeros_matrix<double>(number_of_inequality_constraints_, 1);
                    }
                }

                if (minor_iteration > 1) {
                    temporary_gradient = gradient - temporary_gradient;
                    temporary_parameter = parameter - temporary_parameter;

                    double sc[2];
                    sc[0] = dlib::trans(temporary_parameter) * hessian * temporary_parameter;
                    sc[1] = dlib::trans(temporary_parameter) * temporary_gradient;
                    if (sc[0] * sc[1] > 0) {
                        temporary_parameter = hessian * temporary_parameter;
                        hessian = hessian - (temporary_parameter * dlib::trans(temporary_parameter)) / sc[0] +
                                  (temporary_gradient * dlib::trans(temporary_gradient)) / sc[1];

                    }
                }

                dx = 0.01 * dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_, 1);

                if (lagrangian_parameters_bounded_.second) {
                    gap = dlib::join_rows(
                            dlib::rowm(parameter, dlib::range(0, mm - 1)) -
                            dlib::colm(parameter_bounds, 0),
                            dlib::colm(parameter_bounds, 1) -
                            dlib::rowm(parameter, dlib::range(0, mm - 1))
                    );

                    gap = left_vector_min_right_vector_max(gap);

                    gap = dlib::colm(gap, 0) +
                          std::sqrt(std::numeric_limits<double>::epsilon()) * dlib::ones_matrix<double>(mm, 1);

                    dlib::set_rowm(dx, dlib::range(0, mm - 1)) = pointwise_divide(dlib::ones_matrix<double>(mm, 1),
                                                                                  gap);

                    if (!lagrangian_parameters_bounded_.first) {
                        dlib::set_rowm(dx, dlib::range(mm, number_of_parameters_and_inequality_constraints_ - 1)) =
                                dlib::min(dlib::join_cols(dlib::rowm(dx, dlib::range(0, mm - 1)), dlib::mat(0.01))) *
                                dlib::ones_matrix<double>(number_of_parameters_and_inequality_constraints_ - mm, 1);

                    }
                }
                go = -1.0;
                mu = mu / 10.0;
                dlib::matrix<double> u;
                dlib::matrix<double> cholesky;
                while (go <= 0) {

                    cholesky = dlib::trans(dlib::chol(
                            dlib::make_symmetric(hessian) + mu * dlib::diagm(dlib::pointwise_multiply(dx, dx))));

                    cholesky = dlib::inv_upper_triangular(cholesky);

                    temporary_gradient = dlib::trans(cholesky) * gradient;

                    if (number_of_total_constraints_ == 0) {
                        u = -1 * cholesky * temporary_gradient;

                    } else {
                        // We solve the equation system using QR factorization
                        dlib::qr_decomposition<dlib::matrix<double>> qr(trans(cholesky) * dlib::trans(a));
                        lagrangian_multipliers = qr.solve(temporary_gradient);

                        u = -1 * cholesky *
                            (temporary_gradient - (dlib::trans(cholesky) * dlib::trans(a)) * lagrangian_multipliers);
                    }

                    parameter0 = dlib::rowm(u, dlib::range(0, number_of_parameters_and_inequality_constraints_ - 1)) +
                                 parameter;

                    if (!lagrangian_parameters_bounded_.second) {
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


                if (number_of_inequality_constraints_ != 0) {
                    dlib::set_rowm(cost_vector3,
                                   dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                            dlib::rowm(cost_vector3,
                                       dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) -
                            dlib::subm(pt, dlib::rectangle(2, 0, 2, number_of_inequality_constraints_ - 1));
                }
                if (number_of_total_constraints_ != 0) {
                    dlib::set_rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) =
                            dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) -
                            a * dlib::colm(pt, 2) + b;
                    sob(2) = cost_vector3(0) - dlib::trans(lagrangian_multipliers0) *
                                               dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_)) +
                             rho * std::pow<double>(euclidean_norm(
                                     dlib::rowm(cost_vector3, dlib::range(1, number_of_total_constraints_))), 2);
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

                    if (number_of_inequality_constraints_ != 0) {
                        dlib::set_rowm(cost_vector2,
                                       dlib::range(number_of_equality_constraints_ + 1, number_of_total_constraints_)) =
                                dlib::rowm(cost_vector2, dlib::range(number_of_equality_constraints_ + 1,
                                                                     number_of_total_constraints_)) -
                                dlib::subm(pt, dlib::rectangle(1, 0, 1, number_of_inequality_constraints_ - 1));
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
                    }
                    obm = dlib::max(sob);
                    if (obm < object_function_value) {
                        obn = dlib::min(sob);
                        go = tolerance * (obm - obn) / (object_function_value - obm);
                    }

                    if (sob(1) >= sob(0)) {
                        sob(2) = sob(1);
                        cost_vector3 = cost_vector2;
                        alpha_(2) = alpha_(1);
                        dlib::set_colm(pt, 2) = dlib::colm(pt, 1);
                    } else if (sob(0) <= sob(2)) {
                        sob(2) = sob(1);
                        cost_vector3 = cost_vector2;
                        alpha_(2) = alpha_(1);
                        dlib::set_colm(pt, 2) = dlib::colm(pt, 1);
                    } else {
                        sob(0) = sob(1);
                        cost_vector1 = cost_vector2;
                        alpha_(0) = alpha_(1);
                        dlib::set_colm(pt, 0) = dlib::colm(pt, 1);

                    }

                    if (go >= tolerance) {
                        go = alpha_(2) - alpha_(0);
                    }

                }

                temporary_parameter = parameter;
                temporary_gradient = gradient;
                positive_change = true;
                obn = dlib::min(sob);

                if (object_function_value <= obn) {
                    max_iterations = minor_iteration;
                }

                reduction = (object_function_value - obn) / (1 + std::abs(object_function_value));

                //Reduction too low? Then we end the loop.
                if (reduction < tolerance) {
                    max_iterations = minor_iteration;
                }

                /* Row 262 */
                if (sob(0) < sob(1)) {
                    object_function_value = sob(0);
                    parameter = dlib::colm(pt, 0);
                    cost_vector = cost_vector1;
                } else if (sob(2) < sob(1)) {
                    object_function_value = sob(2);
                    parameter = dlib::colm(pt, 2);
                    cost_vector = cost_vector3;
                } else {
                    object_function_value = sob(1);
                    parameter = dlib::colm(pt, 1);
                    cost_vector = cost_vector2;
                }
            }


            /*
             * Unscale the parameter vector
             * */
            parameter = dlib::pointwise_multiply(
                    parameter,
                    dlib::rowm(scale, dlib::range(number_of_equality_constraints_ + 1,
                                                  number_of_total_constraints_ + number_of_parameters_))
            );


            if (number_of_total_constraints_ != 0) {
                /*
                 * Unscale the lagrange multipliers
                 * */
                lagrangian_multipliers = scale(0) * pointwise_divide(lagrangian_multipliers,
                                                                     dlib::rowm(scale, dlib::range(1,
                                                                                                   number_of_total_constraints_)));

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

            // Guarantee that the hessian matrix is symmetric (in theory, it should already be)
            hessian = dlib::make_symmetric(hessian);

            if (event_log && reduction > tolerance) {
                event_log->push_back(
                        "Warning: Minor optimization routine did not converge. You may need to increase the number of minor iterations.");
            }

        }

    private:
        // Constructor variables
        functor_model &objective_function_;
        const long number_of_parameters_;
        const long number_of_equality_constraints_;
        const long number_of_inequality_constraints_;
        const long number_of_total_constraints_;
        const long number_of_parameters_and_inequality_constraints_;
        const std::pair<bool, bool> &lagrangian_parameters_bounded_;

        // Internal variables
        dlib::matrix<double, 3, 1> alpha_;
    };

} // namespace cppsolnp

#endif //CPP_SOLNP_SUBNP_HPP
