.. _Introduction:

Introduction
============

SOLNP
-----

The SOLNP algorithm by Yinyu Ye (1989) solves the general nonlinear optimization problem below.
In other words, it will find an optimum (possibly local) to the problem below.
The algorithm was originally implemented in `Matlab`_, and have gained some fame through it's R implementation (`RSOLNP`_).
solnp is written in C++ and with the Python wrappers (pysolnp) you have seamless integration with Python, providing high efficiency and ease of use.

.. math::
   :nowrap:

   \begin{gather*}
    \min_{\mathbf{x}} f(\mathbf{x}) \\
    \textrm{s.t.} \\
    \mathbf{g}(\mathbf{x}) = \mathbf{e}_\mathbf{x} \\
    \mathbf{l}_\mathbf{h} < \mathbf{h}(\mathbf{x}) < \mathbf{u}_\mathbf{h} \\
    \mathbf{l}_\mathbf{x} < \mathbf{x} < \mathbf{u}_\mathbf{x} \\
   \end{gather*}

Here:

- :math:`\mathbf{x}` is the optimization parameter.
- :math:`f(\mathbf{x})`,  :math:`\mathbf{h}(\mathbf{x})` and :math:`\mathbf{g}(\mathbf{x})` are smooth constraint functions.
- The constant-valued vector :math:`\mathbf{e}_\mathbf{x}` is the target value for the equality constraint(s).
- The constant-valued vectors :math:`\mathbf{l}_\mathbf{h}` and :math:`\mathbf{u}_\mathbf{h}` are the lower and upper limits resp. for the inequality constraint(s).
- The constant-valued vectors :math:`\mathbf{l}_\mathbf{x}` and :math:`\mathbf{u}_\mathbf{x}` are the lower and upper limits resp. for the parameter constraint(s).

Bold characters indicate the variable or function can be vector-valued. All constraints are optional and vectors can be of any dimension.

.. _RSOLNP: https://cran.r-project.org/web/packages/Rsolnp/index.html
.. _`Matlab`: https://web.stanford.edu/~yyye/matlab/

GOSOLNP
-------

GOSOLNP tries to find the global optimum for the given problem above.
This algorithm was originally implemented in R through `RSOLNP`_.

This is done by:

#. Generate n random starting parameters based on some specified distribution.
#. Evaluate the starting parameters based on one of two evaluation functions (lower value is better):

   a. Objective function :math:`f(\mathbf{x})` for all :math:`\mathbf{x}` that satisfies the inequality constraints :math:`\mathbf{l}_\mathbf{x} < \mathbf{x} < \mathbf{u}_\mathbf{x}`
   b. Penalty Barrier function: :math:`f(\mathbf{x}) + 100 \sum_i(\max(0, 0.9 + l_{x_i} - [\mathbf{g}(\mathbf{x})]_i)^2 + \max(0, 0.9 + [\mathbf{g}(\mathbf{x})]_i - u_{x_i})^2) + \sum_j([h(\mathbf{x})]_j - e_{x_j})^2/100`

#. For the m starting parameters with the lowest evaluation function value, run pysolnp to find nearest optimum.
#. Return the best converged solution among the ones found through the various starting parameters (lowest solution value.)