[![codecov](https://codecov.io/gh/KristerSJakobsson/solnp/branch/master/graph/badge.svg)](https://codecov.io/gh/KristerSJakobsson/solnp)
[![Documentation Status](https://readthedocs.org/projects/solnp/badge/?version=latest)](https://solnp.readthedocs.io/en/latest/?badge=latest)
[![Python Versions](https://img.shields.io/pypi/pyversions/pysolnp.svg)](https://pypi.org/project/pysolnp/)

See full documentation on [http://solnp.readthedocs.io](https://solnp.readthedocs.io/en/latest/).

# pysolnp - Nonlinear optimization with the augmented Lagrange method

## Description
SOLNP solves the general nonlinear optimization problem on the form:
```
    minimize f(x)
      subject to
       g(x) = e_x
   l_h <= h(x) <= u_h
   l_x <   x   < u_X
```
where f(x), g(x) and h(x) are smooth functions.

## Compatibility
Precompiled Wheels are available for CPython:
- Windows: Python 3.6+
- Linux: Python 3.6+
- Mac OS: Python 3.6+

For other systems, or to have BLAS and LAPACK support, please build the wheels manually.
Note: For best results, building it from source is recommended, as BLAS and LAPACK will make a difference.

## Installation
Simply install the package through PyPi with:
`pip install pysolnp`

When compiling from source code you will need CMake.<br>
See the [README for the C++ code](https://github.com/KristerSJakobsson/solnp/blob/master/README.md) for details.

## Usage
Below is the Box example, for the complete example see [/python_examples/example_box.py](/python_examples/example_box.py).
```python
import pysolnp

def f_objective_function(x):
    return -1 * x[0] * x[1] * x[2]

def g_equality_constraint_function(x):
    return [4 * x[0] * x[1] + 2 * x[1] * x[2] + 2 * x[2] * x[0]]

x_starting_point = [1.1, 1.1, 9.0]
x_l = [1.0, 1.0, 1.0]
x_u = [10.0, 10.0, 10.0]
e_x = [100]

result = pysolnp.solve(
    obj_func=f_objective_function,
    par_start_value=x_starting_point,
    par_lower_limit=x_l,
    par_upper_limit=x_u,
    eq_func=g_equality_constraint_function,
    eq_values=e_x)

result.solve_value
result.optimum
result.callbacks
result.converged
```

Output:
```
>>> result.solve_value
-48.11252206814995
>>> result.optimum
[2.8867750707815447, 2.8867750713194273, 5.773407748939196]
>>> result.callbacks
118
>>> result.converged
True
```

## Parameters
The basic signature is:
```python
solve(obj_func: function, par_start_value: List, par_lower_limit: object = None, par_upper_limit: object = None, eq_func: object = None, eq_values: object = None, ineq_func: object = None, ineq_lower_bounds: object = None, ineq_upper_bounds: object = None, rho: float = 1.0, max_major_iter: int = 10, max_minor_iter: int = 10, delta: float = 1e-05, tolerance: float = 0.0001, debug: bool = False) -> pysolnp.Result
```

Inputs:

| Parameter          | Type                      | Default value*   | Description                                                                       |
| -------------------|:--------------------------|:-----------------|-----------------------------------------------------------------------------------|
| obj_func           | Callable\[List, float\]   | -                | The objective function f(x) to minimize.                                          |
| par_start_value    | List                      | -                | The starting parameter x_0.                                                       |
| par_lower_limit    | List                      | None             | The parameter lower limit x_l.                                                    |
| par_upper_limit    | List                      | None             | The parameter upper limit x_u.                                                    |
| eq_func            | Callable\[List, float\]   | None             | The equality constraint function h(x).                                            |
| eq_values          | List                      | None             | The equality constraint values e_x.                                               |
| ineq_func          | Callable\[List, float\]   | None             | The inequality constraint function g(x).                                          |
| ineq_lower_bounds  | List                      | None             | The inequality constraint lower limit g_l.                                        |
| ineq_upper_bounds  | List                      | None             | The inequality constraint upper limit g_l.                                        |
| rho                | float                     | 1.0              | Penalty weighting scalar for infeasability in the augmented objective function.** |
| max_major_iter     | int                       | 400              | Maximum number of outer iterations.                                               |
| max_minor_iter     | int                       | 800              | Maximum number of inner iterations.                                               |
| delta              | float                     | 1e-07            | Step-size for forward differentiation.                                            |
| tolerance          | float                     | 1e-08            | Relative tolerance on optimality.                                                 |
| debug              | bool                      | False            | If set to true some debug output will be printed.                                 |

*Defaults for configuration parameters are based on the defaults for Rsolnp.<br>
**Higher values means the solution will bring the solution into the feasible region with higher weight. Very high values might lead to numerical ill conditioning or slow down convergence.

Output:
The function returns the `pysolnp.Result` with the below properties.

| Property           | Type                  | Description                                           |
| -------------------|:----------------------|-------------------------------------------------------|
| solve_value        | float                 | The value of the objective function at optimum f(x*). |
| optimum            | List\[float\]         | A list of parameters for the optimum x*.              |
| callbacks          | int                   | Number of callbacks done to find this optimum.        |
| converged          | boolean               | Indicates if the algorithm converged or not.          |
| hessian_matrix     | List\[List\[float\]\] | The final Hessian Matrix used by pysolnp.             |

## Use-cases and Applications
* NMPC - Nonlinear model predictive controls-case studies using Matlab, REXYGEN and pysolnp NLP solver under Python environment by Štěpán Ožana. 
[[NMPC Overhead Crane (PDF)](https://github.com/StepanOzana/NMPC/raw/main/NMPC_Overhead_Crane/NMPC_overhead_crane_description.pdf)] 
[[GitHub Source Code](https://github.com/StepanOzana/NMPC)]

## Authors

* **Krister S Jakobsson** - *Implementation* - krister.s.jakobsson@gmail.com

## License

This project is licensed under the Boost License - see the [license](LICENSE.md) file for details.

## Acknowledgments

* **Yinyu Ye** -  Publisher and mastermind behind the original SOLNP algorithm,
[Original Sources](https://web.stanford.edu/~yyye/matlab/)
* **Alexios Ghalanos and Stefan Theussl** - The people behind RSOLNP,
[Github repository](https://github.com/cran/Rsolnp)
* **Davis King** - The mastermind behind Dlib, check out his blog! [Blog](http://blog.dlib.net/) 
