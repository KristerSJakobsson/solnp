//
// Created by krister on 2019/01/20.
//

#ifndef CPP_SOLNP_SOLVE_HPP
#define CPP_SOLNP_SOLVE_HPP


#include "stdafx.h"

namespace cppsolnp {

    namespace {
        static const double sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());

        bool isSquareMatrix(const dlib::matrix<double>& A) {
            return A.nr() == A.nc();
        }
    }


    // Inspired by Matlabs mldivide function
    dlib::matrix<double, 0, 1> solve(const dlib::matrix<double> &A,
                                     dlib::matrix<double, 0, 1> b) {
        double eps = dlib::max(dlib::abs(A)) * sqrtEps / 100;

        if (!isSquareMatrix(A)) {
            // Use QR Solver for non-square matrix A
            return dlib::pinv(A)*b;
        } else {

            bool upper_triangular = true, lower_triangular = true, symmetric = true;

            dlib::matrix<double, 1, 0> row_right_of_diagonal;
            dlib::matrix<double, 0, 1> col_below_diagonal;

            for (int i = 1; i < A.nr() - 1; ++i) {
                // l t r b
                col_below_diagonal = dlib::subm(A, dlib::rectangle(i-1, i, i-1, A.nr()-1));
                row_right_of_diagonal = dlib::subm(A, dlib::rectangle(i, i-1, A.nc()-1, i-1));

                if (dlib::abs(col_below_diagonal) > eps)
                    upper_triangular = false;
                if (dlib::abs(row_right_of_diagonal) > eps)
                    lower_triangular = false;
                if (dlib::abs(col_below_diagonal - dlib::trans(row_right_of_diagonal)) > eps)
                    symmetric = false;

                if (!upper_triangular &&
                    !lower_triangular &&
                    !symmetric)
                    break;
            }

            if (symmetric) {
                dlib::cholesky_decomposition<dlib::matrix<double>> decomposition(A);
                return decomposition.solve(b);
            }

            if (upper_triangular) {
                return dlib::inv_upper_triangular(A) * b;
            }

            if (lower_triangular) {
                return dlib::inv_lower_triangular(A) * b;
            }

            dlib::lu_decomposition<dlib::matrix<double>> decomposition(A);
            return decomposition.solve(b);

        }
    }

}

#endif //CPP_SOLNP_SOLVE_HPP
