"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ./pysolnp" from the cpp_solnp folder.
2) Run this file with Python
"""

import pysolnp


# Alkyla function from the original documentation of SOLNP
# This problem has both equality and inequality constraints, as well as the parameter bounds.


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
        (-0.222 * x[0] + 35.82) / x[8],
        (3.0 * x[6] - 133.0) / x[0],
    ]
    return result


start_point = []
equality_values = [0, 0, 0]
inequality_lower_bounds = [0, 0, 0, 0]
inequality_upper_bounds = [100, 100, 100, 100]

if __name__ == "__main__":
    val_optim, x_optim, debug = pysolnp.solve(x0=start_point,
                                              obj_func=alkyla_objective_function,
                                              eq_func=alkyla_equality_function,
                                              eq_values=equality_values,
                                              ineq_func=alkyla_inequality_function,
                                              ineq_lower_bounds=inequality_lower_bounds,
                                              ineq_upper_bounds=inequality_upper_bounds,
                                              rho=1.0,
                                              max_major_iter=10,
                                              max_minor_iter=10,
                                              delta=1e-5,
                                              tolerance=1e-4,
                                              debug=True)
    result = "Optimal value %d found at %d" % (val_optim, x_optim)
    print(result)
    print(debug)