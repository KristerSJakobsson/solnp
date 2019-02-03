//
// Created by krister on 2019/01/20.
//

#ifndef CPP_SOLNP_UTILS_HPP
#define CPP_SOLNP_UTILS_HPP

#include "stdafx.h"

namespace cppsolnp {

    typedef std::vector<std::string> log_list;
    typedef std::unique_ptr<std::vector<std::string>> log_list_ptr;

    /* Calculates the 2-norm conditional number using SVD.
       The 2-norm conditional number is defined as
        cons2(M) = euclidean_norm(M) * euclidean_norm(M^-1) = sigma_1/sigma_n
        where sigma are the singular values of M, where
        sigma_1 is the biggest singular value and sigma_n is the smallest.
    */
    template<typename M>
    double conditional_number(const M &matrix) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<M>::value);

        /* TODO: This function can be optimized.*/

        int dim_min = std::min(matrix.nc(), matrix.nr());
        int dim_max = std::max(matrix.nc(), matrix.nr());

        dlib::matrix<double> singular_value_matrix;

        long result;
        dlib::matrix<double> dummy_max(dim_max, dim_max);
        dlib::matrix<double> dummy_min(dim_min, dim_min);

        if (matrix.nr() >= matrix.nc()) {
            result = dlib::svd2(false, false, matrix, dummy_max, singular_value_matrix, dummy_min);
        } else {
            result = dlib::svd2(false, false, dlib::trans(matrix), dummy_max, singular_value_matrix, dummy_min);
        }
        if (result != 0) {
            throw std::runtime_error("Error: Singular value decomposition failed.");
            // Alternative: Use the "svd" algorithm.

        }

        // Notably, the svd2 routine will not return the singular values in order!
        double singular_value_min, singular_value_max;
        dlib::find_min_and_max(singular_value_matrix, singular_value_min, singular_value_max);

        return singular_value_max / singular_value_min;
    }

    template<typename V>
    double inline euclidean_norm(V vector) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<V>::value);
        return std::sqrt(dlib::dot(vector, vector));
    }

    template<typename V>
    double inline infinity_norm(V vector) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<V>::value);
        return dlib::max(dlib::abs(vector));
    }

    template<typename M1,
            typename M2>
    dlib::matrix<double> pointwise_divide(const M1 &numerator_matrix, const M2 &denominator_matrix) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<M1>::value);
        COMPILE_TIME_ASSERT(dlib::is_matrix<M2>::value);
        if (numerator_matrix.nc() != denominator_matrix.nc() || numerator_matrix.nr() != denominator_matrix.nr()) {
            throw std::runtime_error("Error: Tried to divide two matrixes of different size.");
        }
        int max_col = numerator_matrix.nc(), max_row = numerator_matrix.nr();
        dlib::matrix<double> result(max_row, max_col);
        double denom;
        for (int row = 0; row < max_row; ++row) {
            for (int col = 0; col < max_col; ++col) {

                if ((denom = denominator_matrix(row, col)) == 0.0) {
                    result(row, col) = std::numeric_limits<double>::infinity();
                } else if (denom == -0.0) {
                    result(row, col) = -std::numeric_limits<double>::infinity();
                } else {
                    result(row, col) = numerator_matrix(row, col) / denom;
                }
            }
        }

        return result;
    }

    /* Max between each element in a vector and a scalar.*/
    template<typename V>
    V elementwise_max(V vector, const double &scalar) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<V>::value);
        COMPILE_TIME_ASSERT(V::NC <= 1);
        if (vector.nc() > 1 || vector.nr() == 0) {
            throw std::runtime_error("Error: Please input a column vector.");
        }

        for (auto element : vector) {
            element = std::max(element, scalar);
        }

        return vector;
    }

    /* Min between each element in a vector and a scalar.*/
    template<typename V>
    V elementwise_min(V vector, const double &scalar) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<V>::value);
        COMPILE_TIME_ASSERT(V::NC <= 1);
        if (vector.nc() > 1 || vector.nr() == 0) {
            throw std::runtime_error("Error: Please input a column vector.");
        }

        for (auto element : vector) {
            element = std::min(element, scalar);
        }
        return vector;
    }

    /* Arranges two column vectors, inputed as a matrix,
    so that the left vector contains the min and
    the right vector contains the max.*/
    template<typename M>
    void two_vector_sort(M &matrix) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<M>::value);
        COMPILE_TIME_ASSERT(M::NC == 2 || M::NC == 0);
        if (matrix.nc() != 2 || matrix.nr() == 0) {
            throw std::runtime_error("Error: Invalid input.");
        }

        dlib::matrix<double, 1, 2> temp;
        for (int row = 0; row < matrix.nr(); ++row) {
            temp(0) = dlib::min(dlib::rowm(matrix, row));
            temp(1) = dlib::max(dlib::rowm(matrix, row));
            dlib::set_rowm(matrix, row) = temp;
        }
        return;
    }

    /* Returns the vector with the max of each row in a matrix.*/
    template<typename M>
    dlib::matrix<double, 0, 1> rowwise_max(M &matrix) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<M>::value);

        if (matrix.nr() == 0 || matrix.nc() == 0) {
            throw std::runtime_error("Error: Invalid input.");
        }

        dlib::matrix<double, 0, 1> return_vector(matrix.nr());
        for (int row = 0; row < matrix.nr(); ++row) {
            return_vector(row) = dlib::max(dlib::rowm(matrix, row));
        }
        return return_vector;
    }


    /* Saves a matrix to a string. If flatten is true, no row breaks will be included.*/
    template<typename M>
    std::string to_string(M matrix, bool flatten = false) {
        COMPILE_TIME_ASSERT(dlib::is_matrix<M>::value);
        std::string return_value = "[";
        for (int row = 0; row < matrix.nr(); ++row) {
            for (int col = 0; col < matrix.nc(); ++col) {
                return_value += std::to_string(matrix(row, col)) + " ";
                if (col < matrix.nc() - 1 && row < matrix.nr() - 1) {
                    return_value += " ";
                } else if (col == matrix.nc() - 1 && row < matrix.nr() - 1) {
                    return_value += ";";
                }
            }
            if (row < matrix.nr() - 1 && flatten == false) {
                return_value += "\n";
            } else if (row < matrix.nr() - 1 && flatten == true) {
                return_value += " ";
            } else {
                return_value += "]";
            }
        }

        return return_value;
    }


}

#endif //CPP_SOLNP_UTILS_HPP
