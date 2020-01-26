# pysolnp - Nonlinear optimization with the augmented Lagrange method

## Description
SOLNP solves the general nonlinear optimization problem on the form:
```
    minimize f(x)
      subject to
       g(x) = e_x
   l_h <= h(x) <= u_x
   l_x <=  x   <= u_X
```
where f(x), g(x) and h(x) are smooth functions.

## Installation

Simply install the package:
`pip install pysolnp`

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
```

Output:
```
>>> result.solve_value
-48.11252206814995
>>> result.optimum
[2.8867750707815447, 2.8867750713194273, 5.773407748939196]
>>> result.callbacks
118
```

## Parameters
The basic signature is:
```python
solve(obj_func: function, par_start_value: list, par_lower_limit: object = None, par_upper_limit: object = None, eq_func: object = None, eq_values: object = None, ineq_func: object = None, ineq_lower_bounds: object = None, ineq_upper_bounds: object = None, rho: float = 1.0, max_major_iter: int = 10, max_minor_iter: int = 10, delta: float = 1e-05, tolerance: float = 0.0001, debug: bool = False) -> pysolnp.Result
```

Inputs:

| Parameter          | Type                      | Default value*   | Description                                                                       |
| -------------------|:--------------------------|:-----------------|-----------------------------------------------------------------------------------|
| obj_func           | Callable\[list, float\]   | -                | The objective function f(x) to minimize.                                          |
| par_start_value    | list                      | -                | The starting parameter x_0.                                                       |
| par_lower_limit    | list                      | None             | The parameter lower limit x_l.                                                    |
| par_upper_limit    | list                      | None             | The parameter upper limit x_u.                                                    |
| eq_func            | Callable\[list, float\]   | None             | The equality constraint function h(x).                                            |
| eq_values          | list                      | None             | The equality constraint values e_x.                                               |
| ineq_func          | Callable\[list, float\]   | None             | The inequality constraint function g(x).                                          |
| ineq_lower_bounds  | list                      | None             | The inequality constraint lower limit g_l.                                        |
| ineq_upper_bounds  | list                      | None             | The inequality constraint upper limit g_l.                                        |
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

| Property           | Type           | Description                                           |
| -------------------|:---------------|-------------------------------------------------------|
| solve_value        | double         | The value of the objective function at optimum f(x*). |
| optimum            | list\[double\] | A list of parameters for the optimum x*.              |
| callbacks          | int            | Number of callbacks done to find this optimum.        |

## Authors

* **Krister S Jakobsson** - *Implementation* - krister.s.jakobsson@gmail.com

## License

This project is licensed under the Boost License - see the [license](LICENSE.md) file for details.

## Acknowledgments

* **Yinyu Ye** -  Publisher of the original algorithm,
[Original Sources](https://web.stanford.edu/~yyye/matlab/)
* **Rsolnp** - An R implementation of SOLNP,
[Github repository](https://github.com/cran/Rsolnp)