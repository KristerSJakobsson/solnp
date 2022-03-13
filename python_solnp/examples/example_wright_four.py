"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ." from the root folder.
2) Run this file with Python
"""

import pysolnp
from math import sqrt


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


def solve_wright_four():
    result = pysolnp.solve(
        obj_func=wright_four_objective_function,
        par_start_value=starting_point,
        eq_func=wright_four_equality_function,
        eq_values=equality_values)
    return result


if __name__ == "__main__":
    result = solve_wright_four()

    final_parameters = result.optimum
    print(final_parameters)
    print(result.solve_value)
    print(result.callbacks)

    print(wright_four_equality_function(final_parameters))

    final_objective_value = wright_four_equality_function(final_parameters)

    equality_constaints = wright_four_equality_function(final_parameters)

    for index, value in enumerate(equality_constaints):
        distance_equality_constraint = value - equality_values[index]
        print("Distance for Equality Constraint index %s: distance %s" % (index, distance_equality_constraint))
