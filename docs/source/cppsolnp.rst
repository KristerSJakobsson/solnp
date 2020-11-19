.. _C++ solnp:

C++ solnp
=========
C++ solnp implements the SOLNP algorithm explained in the :ref:`Introduction<Introduction>` section.

Installation
------------
C++ solnp is a header-only library, so the only thing you need to do is:

- Add dlib to your project.
- Include the C++ solnp header :code:`solnp.hpp` in your project.

Note that the CMake script in the root of the repository is mainly intended to compile wrappers and unit tests.

Methods
-------
The solnp function can be called both with and without Hessian Matrix, as shown below.

.. code-block:: c

    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    double cppsolnp::solnp(
                functor_model functor,
                parameter_input &parameter_data,
                const inequality_constraint_vectors &inequality_constraint_data,
                const std::shared_ptr<std::vector<std::string>> &event_log = nullptr,
                double rho = 1.0,
                int maximum_major_iterations = 400,
                int maximum_minor_iterations = 800,
                const double &delta = 1e-7,
                const double &tolerance = 1e-8
        )

.. code-block:: c
   :emphasize-lines: 9

    template<
            typename functor_model,
            typename parameter_input,
            typename inequality_constraint_vectors>
    double cppsolnp::solnp(
                functor_model functor,
                parameter_input &parameter_data,
                const inequality_constraint_vectors &inequality_constraint_data,
                dlib::matrix<double> &hessian_matrix,
                const std::shared_ptr<std::vector<std::string>> &event_log = nullptr,
                double rho = 1.0,
                int maximum_major_iterations = 400,
                int maximum_minor_iterations = 800,
                const double &delta = 1e-7,
                const double &tolerance = 1e-8
        )

Templates:

- functor_model: Represents a callable functor or lambda.
- parameter_input: Represents the parameter input, this should be a dynamic-height or static-height dlib::matrix with:

   - 1 column representing the start point :math:`\mathbf{x*}` if the problem has no parameter bounds.
   - 2 columns representing the lower limit :math:`\mathbf{x_l}` and upper limit :math:`\mathbf{x_u}`, in this case th inital guess will be the middle point of the bounds.
   - 3 columns representing the start point :math:`\mathbf{x*}`, lower limit :math:`\mathbf{x_l}` and upper limit :math:`\mathbf{x_u}`.

- inequality_constraint_vectors: Represents the inequality constraint, this should be a dynamic-height or static-height dlib::matrix with:

   - 2 columns representing the lower bounds :math:`\mathbf{l}_\mathbf{x}` and upper bounds :math:`\mathbf{u}_\mathbf{x}`.

Note: In contrast to pysolnp, C++ solnp will assume that the equality constraints equal 0, so supply a function :math:`\mathbf{g}_{new}(\mathbf{x}) = \mathbf{g}(\mathbf{x}) - \mathbf{e}_\mathbf{x} = \mathbf{0} \\`.

Results
-------
The solnp function will return the solve result value as a double.
The parameter_data matrix will be modified in-memory and contain the corresponding solution.
If a logger is supplied as a smart-pointer to an std::vector of strings, various log messages will be stored in it.

Example 1: Box Problem
----------
The Box Problem is a common example function used for testing optimization algorithms.
It has one equality constraint and variable bounds.

The below code is taken from the C++ solnp unit tests.

.. code-block:: c

    #include "catch.hpp"
    #include "solnp.hpp"

    dlib::matrix<double, 2, 1> box(const dlib::matrix<double, 3, 1> &m)
    {
        const double x1 = m(0);
        const double x2 = m(1);
        const double x3 = m(2);

        dlib::matrix<double, 2, 1> return_values(2);
        // Function value
        return_values(0) = -1 * x1 * x2 * x3;
        // Equality constraint
        return_values(1) = 4 * x1 * x2 + 2 * x2 * x3 + 2 * x3 * x1 - 100;
        return return_values;
    }


    struct box_functor {
    public:
        box_functor() = default;;

        dlib::matrix<double, 2, 1> operator()(const dlib::matrix<double, 3, 1> &x) {
            return box(x);
        }
    };

    TEST_CASE("Optimize the Box function", "[box]") {

        dlib::matrix<double, 3, 3> parameter_data;
        parameter_data =
                1.1, 1.0, 10.0,
                1.1, 1.0, 10.0,
                9.0, 1.0, 10.0;

        dlib::matrix<double, 0, 0> ib;

        std::shared_ptr<std::vector<std::string>> logger = std::make_shared<std::vector<std::string>>();

        double calculate = cppsolnp::solnp(box_functor(), parameter_data, ib, logger, 1.0, 10, 10, 1e-5, 1e-4);

        dlib::matrix<double, 0, 1> result = dlib::colm(parameter_data, 0);

        // Check the parameters
        CHECK(result(0) == Approx(2.886775069536727));
        CHECK(result(1) == Approx(2.886775072009683));
        CHECK(result(2) == Approx(5.773407750048355));

        REQUIRE(calculate <= -48.112522068150462);

    }
