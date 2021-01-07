.. _Introduction:

Introduction
=============

The SOLNP algorithm by Yinyu Ye (1989) solves the general nonlinear optimization problem below.
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