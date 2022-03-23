"""
To test this algorithm, then:
1) Install this package, for example through pip or by running "pip install ." from the root folder.
2) Run this file with Python
"""

import pysolnp


simple_quadratic = lambda x: x[0] * x[0]

starting_points = [1]

def solve_simple_quadratic():
    result = pysolnp.solve(
        obj_func=simple_quadratic,
        par_start_value=starting_points)
    return result


if __name__ == "__main__":
    result = solve_simple_quadratic()

    final_parameters = result.optimum
    print(final_parameters)
    print(result.solve_value)
    print(result.callbacks)
    print(result.converged)
