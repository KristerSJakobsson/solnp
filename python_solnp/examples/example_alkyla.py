"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ." from the root folder.
2) Run this file with Python
"""

import pysolnp


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


def solve_alkyla():
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
        rho=0.0)
    return result


if __name__ == "__main__":
    result = solve_alkyla()

    final_parameters = result.optimum
    print(final_parameters)
    print(result.solve_value)
    print(result.callbacks)

    print(alkyla_equality_function(final_parameters))
    print(alkyla_inequality_function(final_parameters))

    final_objective_value = alkyla_objective_function(final_parameters)

    equality_constaints = alkyla_equality_function(final_parameters)

    inequality_constraints = alkyla_inequality_function(final_parameters)

    for index, value in enumerate(final_parameters):
        distance_to_lower = value - lower_bound[index]
        distance_to_over = upper_bound[index] - value
        print("Distance for Parameter Constraint index %s: lower %s upper %s" % (
        index, distance_to_lower, distance_to_over))

    for index, value in enumerate(inequality_constraints):
        distance_to_lower = value - inequality_lower_bounds[index]
        distance_to_over = inequality_upper_bounds[index] - value
        print("Distance for Inequality Constraint index %s: lower %s upper %s" % (
        index, distance_to_lower, distance_to_over))

    for index, value in enumerate(equality_constaints):
        distance_equality_constraint = value - equality_values[index]
        print("Distance for Equality Constraint index %s: distance %s" % (index, distance_equality_constraint))
