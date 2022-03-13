"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ." from the root folder.
2) Run this file with Python
"""

import pysolnp


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


def solve_box():
    result = pysolnp.solve(
        obj_func=box_objective_function,
        par_start_value=starting_point,
        par_lower_limit=lower_bound,
        par_upper_limit=upper_bound,
        eq_func=box_equality_function,
        eq_values=equality_values)
    return result


if __name__ == "__main__":
    result = solve_box()

    final_parameters = result.optimum
    print(final_parameters)
    print(result.solve_value)
    print(result.callbacks)

    print(box_equality_function(final_parameters))

    final_objective_value = box_objective_function(final_parameters)

    equality_constaints = box_equality_function(final_parameters)

    for index, value in enumerate(final_parameters):
        distance_to_lower = value - lower_bound[index]
        distance_to_over = upper_bound[index] - value
        print("Distance for Parameter Constraint index %s: lower %s upper %s" % (
        index, distance_to_lower, distance_to_over))

    for index, value in enumerate(equality_constaints):
        distance_equality_constraint = value - equality_values[index]
        print("Distance for Equality Constraint index %s: distance %s" % (index, distance_equality_constraint))
