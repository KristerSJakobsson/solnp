import unittest
import pysolnp
from math import sqrt


class TestExtension(unittest.TestCase):

    def test_solve_simple_quadratic(self):
        simple_quadratic = lambda x: x[0] * x[0]

        starting_points = [1]

        result = pysolnp.solve(
            obj_func=simple_quadratic,
            par_start_value=starting_points)

        self.assertTrue(result.converged)

        # y=x^2 has minimum for x=0, y=0
        self.assertAlmostEqual(result.optimum[0], 0.0)
        self.assertAlmostEqual(result.solve_value, 0.0)

    def test_fail_impossible_problem_with_inequality_constraints(self):
        # This represents the impossible problem:
        # min f(x) = x^2
        #   subject to
        #   0 < x < 1
        #   3 < x < 4
        # Notably, the constraints imply there is no solution.

        simple_quadratic = lambda x: x[0] * x[0]

        starting_points = [1]

        def inequality_constraints(x):
            result = [
                x[0],
                x[0]
            ]
            return result

        result = pysolnp.solve(
            obj_func=simple_quadratic,
            par_start_value=starting_points,
            ineq_func=inequality_constraints,
            ineq_lower_bounds=[0, 3],
            ineq_upper_bounds=[1, 4])

        self.assertFalse(result.converged)

    def test_solve_alkyla(self):
        def alkyla_objective_function(x):
            result = -0.63 * x[3] * x[6] + 50.4 * x[0] + 3.5 * x[1] + x[2] + 33.6 * x[4]
            return result

        def alkyla_equality_function(x):
            result = [
                98.0 * x[2] - 0.1 * x[3] * x[5] * x[8] - x[2] * x[5],
                1000.0 * x[1] + 100.0 * x[4] - 100.0 * x[0] * x[7],
                122.0 * x[3] - 100.0 * x[0] - 100.0 * x[4],
            ]
            return result

        def alkyla_inequality_function(x):
            result = [
                (1.12 * x[0] + 0.13167 * x[0] * x[7] - 0.00667 * x[0] * x[7] * x[7]) / x[3],
                (1.098 * x[7] - 0.038 * x[7] * x[7] + 0.325 * x[5] + 57.25) / x[6],
                (-0.222 * x[9] + 35.82) / x[8],
                (3.0 * x[6] - 133.0) / x[9],
            ]
            return result

        starting_point = [17.45, 12.0, 110.0, 30.0, 19.74, 89.2, 92.8, 8.0, 3.6, 155.0]
        lower_bound = [0.0, 0.0, 0.0, 10.0, 0.0, 85.0, 10.0, 3.0, 1.0, 145.0]
        upper_bound = [20.0, 16.0, 120.0, 50.0, 20.0, 93.0, 95.0, 12.0, 4.0, 162.0]
        equality_values = [0, 0, 0]
        inequality_lower_bounds = [0, 0, 0, 0]
        inequality_upper_bounds = [100, 100, 100, 100]

        result = pysolnp.solve(
            obj_func=alkyla_objective_function,
            par_start_value=starting_point,
            par_lower_limit=lower_bound,
            par_upper_limit=upper_bound,
            eq_func=alkyla_equality_function,
            eq_values=equality_values,
            ineq_func=alkyla_inequality_function,
            ineq_lower_bounds=inequality_lower_bounds,
            ineq_upper_bounds=inequality_upper_bounds,
            rho=0.0,
            max_major_iter=10,
            max_minor_iter=10,
            delta=1e-5,
            tolerance=1e-4)

        final_parameters = result.optimum
        self.assertTrue(result.converged)

        equality_constaints = alkyla_equality_function(final_parameters)

        inequality_constraints = alkyla_inequality_function(final_parameters)

        for index, value in enumerate(final_parameters):
            self.assertLessEqual(value, upper_bound[index] + 1e-5)
            self.assertGreaterEqual(value, lower_bound[index] - 1e-5)

        for index, value in enumerate(inequality_constraints):
            self.assertLessEqual(value, inequality_upper_bounds[index] + 1e-5)
            self.assertGreaterEqual(value, inequality_lower_bounds[index] - 1e-5)

        for index, value in enumerate(equality_constaints):
            self.assertAlmostEqual(value, equality_values[index], 2)

        self.assertLessEqual(result.solve_value,
                             -1.726412486481025e2)  # This reference value is taken form the Matlab reference tests

    def test_solve_box(self):

        def box_objective_function(x):
            result = -1 * x[0] * x[1] * x[2]
            return result

        def box_equality_function(x):
            result = [
                4 * x[0] * x[1] + 2 * x[1] * x[2] + 2 * x[2] * x[0]
            ]
            return result

        starting_point = [1.1,
                          1.1,
                          9.0]
        lower_bound = [1.0,
                       1.0,
                       1.0]
        upper_bound = [10.0,
                       10.0,
                       10.0]

        equality_values = [100]

        result = pysolnp.solve(
            obj_func=box_objective_function,
            par_start_value=starting_point,
            par_lower_limit=lower_bound,
            par_upper_limit=upper_bound,
            eq_func=box_equality_function,
            eq_values=equality_values)

        final_parameters = result.optimum
        self.assertTrue(result.converged)

        equality_constaints = box_equality_function(final_parameters)

        for index, value in enumerate(final_parameters):
            self.assertLessEqual(value, upper_bound[index] + 1e-5)
            self.assertGreaterEqual(value, lower_bound[index] - 1e-5)

        for index, value in enumerate(equality_constaints):
            self.assertAlmostEqual(value, equality_values[index], 2)

        self.assertLessEqual(result.solve_value,
                             -48.112480408240664)  # This reference value is taken form the Matlab reference tests

    def test_solve_rosen_suzuki(self):

        def rosen_suzuki_objective_function(x):
            result = x[0] * x[0] + x[1] * x[1] + 2 * x[2] * x[2] + x[3] * x[3] - 5 * x[0] - 5 * x[1] - 21 * x[2] + 7 * \
                     x[3]
            return result

        def rosen_suzuki_inequality_function(x):
            result = [
                8 - x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3] - x[0] + x[1] - x[2] + x[3],
                10 - x[0] * x[0] - 2 * x[1] * x[1] - x[2] * x[2] - 2 * x[3] * x[3] + x[0] + x[3],
                5 - 2 * x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - 2 * x[0] + x[1] + x[3],
            ]
            return result

        starting_point = [1, 1, 1, 1]
        inequality_lower_bounds = [0, 0, 0]
        inequality_upper_bounds = [1000, 1000, 1000]

        result = pysolnp.solve(
            obj_func=rosen_suzuki_objective_function,
            par_start_value=starting_point,
            ineq_func=rosen_suzuki_inequality_function,
            ineq_lower_bounds=inequality_lower_bounds,
            ineq_upper_bounds=inequality_upper_bounds,
            rho=1.0,
            max_major_iter=10,
            max_minor_iter=10,
            delta=1e-5,
            tolerance=1e-4)

        final_parameters = result.optimum
        self.assertTrue(result.converged)

        inequality_constraints = rosen_suzuki_inequality_function(final_parameters)

        for index, value in enumerate(inequality_constraints):
            self.assertLessEqual(value, inequality_upper_bounds[index] + 1e-5)
            self.assertGreaterEqual(value, inequality_lower_bounds[index] - 1e-5)

        self.assertLessEqual(result.solve_value,
                             -43.999)  # This reference value is taken form the C++ reference tests

    def test_solve_wright_four(self):

        def wright_four_objective_function(x):
            x1 = x[0]
            x2 = x[1]
            x3 = x[2]
            x4 = x[3]
            x5 = x[4]
            result = (x1 - 1) * (x1 - 1) + \
                     (x1 - x2) * (x1 - x2) + \
                     (x2 - x3) * (x2 - x3) * (x2 - x3) + \
                     (x3 - x4) * (x3 - x4) * (x3 - x4) * (x3 - x4) + \
                     (x4 - x5) * (x4 - x5) * (x4 - x5) * (x4 - x5)
            return result

        def wright_four_equality_function(x):
            x1 = x[0]
            x2 = x[1]
            x3 = x[2]
            x4 = x[3]
            x5 = x[4]
            result = [
                x1 + x2 * x2 + x3 * x3 * x3 - 2 - 3 * sqrt(2.0),
                x2 - x3 * x3 + x4 + 2 - 2 * sqrt(2.0),
                x1 * x5 - 2
            ]
            return result

        starting_point = [1.0,
                          1.0,
                          1.0,
                          1.0,
                          1.0]

        equality_values = [0.0,
                           0.0,
                           0.0]

        result = pysolnp.solve(
            obj_func=wright_four_objective_function,
            par_start_value=starting_point,
            eq_func=wright_four_equality_function,
            eq_values=equality_values)

        final_parameters = result.optimum
        self.assertTrue(result.converged)

        equality_constaints = wright_four_equality_function(final_parameters)

        for index, value in enumerate(equality_constaints):
            self.assertAlmostEqual(value, equality_values[index], 2)

        self.assertLessEqual(result.solve_value,
                             0.029310831002758)  # This reference value is taken form the Matlab reference tests


if __name__ == "__main__":
    nose.main()
