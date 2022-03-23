"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ." from the root folder.
2) Run this file with Python
"""

import pysolnp
import math


def powell_objective_function(x):
    result = pow(math.e, x[0] * x[1] * x[2] * x[3] * x[4])
    return result


def powell_equality_function(x):
    z1 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4]
    z2 = x[1] * x[2] - 5 * x[3] * x[4]
    z3 = x[0] * x[0] * x[0] + x[1] * x[1] * x[1]
    return [z1, z2, z3]


starting_point = [-2, 2, 2, -1, -1]
equality_values = [10, 0, -1]


def solve_powell():
    result = pysolnp.solve(
        obj_func=powell_objective_function,
        par_start_value=starting_point)
    # eq_func=powell_equality_function,
    # eq_values=equality_values)
    return result


if __name__ == "__main__":
    result = solve_powell()

    final_parameters = result.optimum
    print(final_parameters)
    print(result.solve_value)
    print(result.callbacks)

    print(powell_equality_function(final_parameters))

    final_objective_value = powell_objective_function(final_parameters)

    equality_constaints = powell_equality_function(final_parameters)

    for index, value in enumerate(equality_constaints):
        distance_equality_constraint = value - equality_values[index]
        print("Distance for Equality Constraint index %s: distance %s" % (index, distance_equality_constraint))
