Python SOLNP (pysolnp)
======================

.. image:: https://codecov.io/gh/KristerSJakobsson/solnp/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/KristerSJakobsson/solnp
   :alt: Codecov Status pysolnp
.. image:: https://img.shields.io/pypi/pyversions/pysolnp.svg
    :target: https://pypi.org/project/pysolnp/

pysolnp provides Python with the power of the SOLNP algorithm explained in :ref:`Introduction<Introduction>` section.
It is simply a Python wrapper for :ref:`C++ solnp<C++ solnp>`.

Installation
------------

In most situations, installing with the package installer for Python, pip, will work:

::

    $ pip install pysolnp

Precompiled Wheels are available for CPython:

- Windows: Python 3.6+
- Linux: Python 3.6+
- Mac OS: Python 3.6+

For other systems, or to have BLAS and LAPACK support, please build the wheels manually.
::

    $ pip install --no-binary :all: pysolnp

Note that this requires CMake.

Method
------

.. code-block:: python

    solve(obj_func: function,
          par_start_value: List,
          par_lower_limit: object = None,
          par_upper_limit: object = None,
          eq_func: object = None,
          eq_values: object = None,
          ineq_func: object = None,
          ineq_lower_bounds: object = None,
          ineq_upper_bounds: object = None,
          rho: float = 1.0,
          max_major_iter: int = 10,
          max_minor_iter: int = 10,
          delta: float = 1e-05,
          tolerance: float = 0.0001,
          debug: bool = False) -> pysolnp.Result

Inputs
-------

+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| Parameter          | Type                      | Default value [#note1]_  | Description                                                                               |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| obj_func           | Callable\[List, float\]   |                          | The objective function :math:`f(x)` to minimize.                                          |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| par_start_value    | List                      |                          | The starting parameter :math:`x_0`.                                                       |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| par_lower_limit    | List                      | None                     | The parameter lower limit :math:`x_l`.                                                    |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| par_upper_limit    | List                      | None                     | The parameter upper limit :math:`x_u`.                                                    |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| eq_func            | Callable\[List, float\]   | None                     | The equality constraint function :math:`h(x)`.                                            |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| eq_values          | List                      | None                     | The equality constraint values :math:`e_x`.                                               |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| ineq_func          | Callable\[List, float\]   | None                     | The inequality constraint function :math:`g(x)`.                                          |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| ineq_lower_bounds  | List                      | None                     | The inequality constraint lower limit :math:`g_l`.                                        |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| ineq_upper_bounds  | List                      | None                     | The inequality constraint upper limit :math:`g_l`.                                        |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| rho                | float                     | 1.0                      | Penalty weighting scalar for infeasability in the augmented objective function. [#note2]_ |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| max_major_iter     | int                       | 400                      | Maximum number of outer iterations.                                                       |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| max_minor_iter     | int                       | 800                      | Maximum number of inner iterations.                                                       |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| delta              | float                     | 1e-07                    | Step-size for forward differentiation.                                                    |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| tolerance          | float                     | 1e-08                    | Relative tolerance on optimality.                                                         |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+
| debug              | bool                      | False                    | If set to true some debug output will be printed.                                         |
+--------------------+---------------------------+--------------------------+-------------------------------------------------------------------------------------------+

.. [#note1] Defaults for configuration parameters are based on the defaults for Rsolnp.
.. [#note2] Higher values means the solution will bring the solution into the feasible region with higher weight. Very high values might lead to numerical ill conditioning or slow down convergence.

Outputs
-------

The function returns the :code:`pysolnp.Result` object with the below properties.

+--------------------+-----------------------+---------------------------------------------------------------+
| Property           | Type                  | Description                                                   |
+--------------------+-----------------------+---------------------------------------------------------------+
| solve_value        | float                 | The value of the objective function at optimum :math:`f(x*)`. |
+--------------------+-----------------------+---------------------------------------------------------------+
| optimum            | List\[float\]         | A list of parameters for the optimum :math:`x*`.              |
+--------------------+-----------------------+---------------------------------------------------------------+
| callbacks          | int                   | Number of callbacks done to find this optimum.                |
+--------------------+-----------------------+---------------------------------------------------------------+
| converged          | boolean               | Indicates if the algorithm converged or not.                  |
+--------------------+-----------------------+---------------------------------------------------------------+
| hessian_matrix     | List\[List\[float\]\] | The final Hessian Matrix used by pysolnp.                     |
+--------------------+-----------------------+---------------------------------------------------------------+

Example 1: Box Problem
------------------------
The Box Problem is a common example function used for testing optimization algorithms.
It has one equality constraint and variable bounds.

.. code-block:: python

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

Running this will yield the output:

::

    >>> result.solve_value
    -48.11252206814995
    >>> result.optimum
    [2.8867750707815447, 2.8867750713194273, 5.773407748939196]
    >>> result.callbacks
    118
    >>> result.converged
    True

Use-cases and Applications
--------------------------
* NMPC - Nonlinear model predictive controls-case studies using Matlab, REXYGEN and pysolnp NLP solver under Python environment by Štěpán Ožana. [`NMPC Overhead Crane (PDF)`_] [`GitHub Source Code`_] [`Štěpán's Homepage`_]

.. _`NMPC Overhead Crane (PDF)`: https://github.com/StepanOzana/NMPC/raw/main/NMPC_Overhead_Crane/NMPC_overhead_crane_description.pdf
.. _`GitHub Source Code`: https://github.com/StepanOzana/NMPC
