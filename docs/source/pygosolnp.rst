Python GOSOLNP (pygosolnp)
==========================
.. image:: https://codecov.io/gh/KristerSJakobsson/pygosolnp/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/KristerSJakobsson/pygosolnp
   :alt: Codecov Status pygosolnp
.. image:: https://img.shields.io/pypi/pyversions/pytest.svg
    :target: https://pypi.org/project/pytest/

pygosolnp provides Python with the power of the GOSOLNP algorithm explained in :ref:`Introduction<Introduction>` section.
It works as an extension on top of pysolnp by solving the problem multiple times from a randomized set of starting points. This library is implemented purely in Python.

Installation
------------
Works on any environment that supports pysolnp and has Python 3.6+ installed.

Method
------

.. code-block:: python

    pygosolnp.solve(
          obj_func: Callable,
          par_lower_limit: List[float],
          par_upper_limit: List[float],
          eq_func: Optional[Callable] = None,
          eq_values: Optional[List[float]] = None,
          ineq_func: Optional[Callable] = None,
          ineq_lower_bounds: Optional[List[float]] = None,
          ineq_upper_bounds: Optional[List[float]] = None,
          number_of_restarts: int = 1,
          number_of_simulations: int = 20000,
          number_of_processes: Optional[int] = None,
          start_guess_sampling: Union[None, List[Distribution], Sampling] = None,
          seed: Union[None, int] = None,
          evaluation_type: Union[EvaluationType, int] = EvaluationType.OBJECTIVE_FUNC_EXCLUDE_INEQ,
          pysolnp_rho: float = 1.0,
          pysolnp_max_major_iter: int = 10,
          pysolnp_max_minor_iter: int = 10,
          pysolnp_delta: float = 1e-05,
          pysolnp_tolerance: float = 0.0001,
          debug: bool = False) -> Results

Inputs
-------

+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter              | Type                             | Default value [#note1]_                    | Description                                                                                                                                       |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| obj_func               | Callable\[List, float\]          |                                            | The objective function :math:`f(x)` to minimize.                                                                                                  |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| par_lower_limit        | List                             | None                                       | The parameter lower limit :math:`x_l`.                                                                                                            |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| par_upper_limit        | List                             | None                                       | The parameter upper limit :math:`x_u`.                                                                                                            |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| eq_func                | Callable\[List, float\]          | None                                       | The equality constraint function :math:`h(x)`.                                                                                                    |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| eq_values              | List                             | None                                       | The equality constraint values :math:`e_x`.                                                                                                       |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ineq_func              | Callable\[List, float\]          | None                                       | The inequality constraint function :math:`g(x)`.                                                                                                  |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ineq_lower_bounds      | List                             | None                                       | The inequality constraint lower limit :math:`g_l`.                                                                                                |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ineq_upper_bounds      | List                             | None                                       | The inequality constraint upper limit :math:`g_l`.                                                                                                |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| number_of_restarts     | int                              | 1                                          | The `number_of_restarts` best evaluation results are used to run pysolnp `number_of_restarts` times.                                              |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| number_of_simulations  | int                              | 20000                                      | Sets how many randomly generated starting guesses we generate and evaluate with the evaluation function.                                          |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| number_of_processes    | int                              | None                                       | Sets how many parallel processes to run when solving the problem. If None the problem is solved in the main processes.                            |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| start_guess_sampling   | List\[Distribution\] or Sampling | None                                       | A list of distributions for generating starting values, one distribution for each parameter. If None, the Uniform distribution is used. [#note3]_ |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| seed                   | int                              | None                                       | By default the MT19937 Generator is used with timestamp-seed. Optionally an integer seed can be supplied.                                         |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| evaluation_type        | EvaluationType or int            | EvaluationType.OBJECTIVE_FUNC_EXCLUDE_INEQ | Selects the evaluation type from the pygosolnp.EvaluationType enum.                                                                               |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| pysolnp_rho            | float                            | 1.0                                        | pysolnp parameter:Penalty weighting scalar for infeasability in the augmented objective function. [#note2]_                                       |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| pysolnp_max_major_iter | int                              | 400                                        | pysolnp parameter:Maximum number of outer iterations.                                                                                             |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| pysolnp_max_minor_iter | int                              | 800                                        | pysolnp parameter:Maximum number of inner iterations.                                                                                             |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| pysolnp_delta          | float                            | 1e-07                                      | pysolnp parameter:Step-size for forward differentiation.                                                                                          |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| pysolnp_tolerance      | float                            | 1e-08                                      | pysolnp parameter:Relative tolerance on optimality.                                                                                               |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| debug                  | bool                             | False                                      | If set to true some debug output will be printed.                                                                                                 |
+------------------------+----------------------------------+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

.. [#note1] Defaults for configuration parameters are based on the defaults for Rsolnp.
.. [#note2] Higher values means the solution will bring the solution into the feasible region with higher weight. Very high values might lead to numerical ill conditioning or slow down convergence.
.. [#note3] Supply an instance of a class that inherits the abstract class `pygosolnp.sampling.Sampling` to provide starting guesses, see `Example 3: Truncated Normal Distribution`_ and `Example 4: Grid Sampling`_.


Outputs
-------

The function returns the :code:`pygosolnp.Results` object with the below properties.

+--------------------+------------------+---------------------------------------------------------------+
| Property           | Type             | Description                                                   |
+--------------------+------------------+---------------------------------------------------------------+
| best_solution      | Optional[Result] | The best local optimum found for the problem.                 |
+--------------------+------------------+---------------------------------------------------------------+
| all_results        | List\[Result\]   | All restarts and their corresponding local optimum.           |
+--------------------+------------------+---------------------------------------------------------------+
| starting_guesses   | List\[float\]    | All the randomized starting parameters.                       |
+--------------------+------------------+---------------------------------------------------------------+

Each named tuple :code:`pygosolnp.Result` has the below properties.

+--------------------+------------------+---------------------------------------------------------------+
| Property           | Type             | Description                                                   |
+--------------------+------------------+---------------------------------------------------------------+
| obj_value          | float            | The value of the objective function at optimum :math:`f(x*)`. |
+--------------------+------------------+---------------------------------------------------------------+
| parameters         | List\[float\]    | A list of parameters for the local optimum :math:`x*`.        |
+--------------------+------------------+---------------------------------------------------------------+
| converged          | bool             | Boolean which indicates if the solution is within bounds.     |
+--------------------+------------------+---------------------------------------------------------------+

Example 1: Electron Optimization Problem
-----------------------------------------

This is a common benchmark problem for Global Optimization that finds the equilibrium state distribution for  electrons positioned on a conducting sphere.
See the `COPS benchmarking suite`_ for details.

.. _`COPS benchmarking suite`: https://www.mcs.anl.gov/~more/cops/

See full source code on GitHub `/python_examples/example_electron.py`_

.. _`/python_examples/example_electron.py`: https://github.com/KristerSJakobsson/pygosolnp/blob/main/python_examples/example_electron.py

.. code-block:: python

    import pygosolnp
    from math import sqrt
    import time

    number_of_charges = 25


    def obj_func(data):
        x = data[0:number_of_charges]
        y = data[number_of_charges:2 * number_of_charges]
        z = data[2 * number_of_charges:3 * number_of_charges]

        result = 0.0
        for i in range(0, number_of_charges - 1):
            for j in range(i + 1, number_of_charges):
                result += 1.0 / sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)

        return result


    def eq_func(data):
        x = data[0:number_of_charges]
        y = data[number_of_charges:2 * number_of_charges]
        z = data[2 * number_of_charges:3 * number_of_charges]
        result = [None] * number_of_charges
        for i in range(0, number_of_charges):
            result[i] = x[i] ** 2 + y[i] ** 2 + z[i] ** 2

        return result


    parameter_lower_bounds = [-1] * number_of_charges * 3
    parameter_upper_bounds = [1] * number_of_charges * 3

    equality_constraints = [1] * number_of_charges

    if __name__ == '__main__':
        start = time.time()

        results = pygosolnp.solve(
            obj_func=obj_func,
            eq_func=eq_func,
            eq_values=equality_constraints,
            par_lower_limit=parameter_lower_bounds,
            par_upper_limit=parameter_upper_bounds,
            number_of_simulations=20000,  # This represents the number of starting guesses to use
            number_of_restarts=20,  # This specifies how many restarts to run from the best starting guesses
            number_of_processes=None,  # None here means to run everything single-processed
            seed=443,  # Seed for reproducibility, if omitted the default random seed is used (typically cpu clock based)
            pysolnp_max_major_iter=100,  # Pysolnp property
            debug=False)

        end = time.time()

        all_results = results.all_results
        print("; ".join([f"Solution {index + 1}: {solution.obj_value}" for index, solution in enumerate(all_results)]))
        best_solution = results.best_solution
        print(f"Best solution {best_solution.obj_value} for parameters {best_solution.parameters}.")
        print(f"Elapsed time: {end - start} s")


::

    Solution 1: 244.1550118432253; Solution 2: 243.9490050190484; Solution 3: 185.78533081425041; Solution 4: 244.07921194485854; Solution 5: 216.19236253370485; Solution 6: 194.1742137471891; Solution 7: 258.6157748268509; Solution 8: 205.72538678938517; Solution 9: 244.0944480181356; Solution 10: 217.4090464122706; Solution 11: 201.58045387715478; Solution 12: 247.70691375326325; Solution 13: 243.92615570955812; Solution 14: 192.3944392661305; Solution 15: 243.93657263760585; Solution 16: 247.17924771908508; Solution 17: 244.06529702108125; Solution 18: 244.29427536763717; Solution 19: 199.69130383979302; Solution 20: 243.99315264179037
    Best solution 243.92615570955812 for parameters [0.8726149386907173, 0.1488320711741995, -0.8215181712229778, 0.8597822831494584, -0.265961670940264, -0.6664127144955102, -0.6029702658967409, 0.2867960203292267, -0.04380531711098636, 0.9519854892760677, -0.39592769694574026, -0.2160514547351913, -0.21416235954836016, 0.4338472533837847, -0.9411378567701716, 0.6418976636970082, 0.014864034847848012, 0.6981416769347426, 0.4413252856284809, -0.5267725521555819, -0.9148568048943023, -0.5831731928212042, 0.47570915153781534, 0.4089885176760918, 0.008471540399374077, -0.36287443863890595, 0.8618964461129363, 0.5476494687199884, -0.3309316231117961, 0.9582851670742292, -0.6505818085537286, 0.2793946112676732, -0.7596998666078645, 0.65142774983249, 0.30572406841664945, -0.1736400992779951, -0.2357569641249718, -0.9762296783338298, 0.8894784482368485, -0.21768032982807542, 0.44966067028074935, 0.359898210796523, 0.3932146838134686, -0.25429503229562933, -0.6621520897149067, 0.0002565729867240561, 0.6081775900274631, -0.8731755460834034, -0.07630776960802095, -0.7462707639808169, 0.32690759610807246, 0.4847543563757037, -0.15870866693945487, -0.38892531575475037, -0.10466177783304143, 0.36421374544164403, -0.7472412325499505, -0.583622807257543, -0.7574487346380878, -0.01614470971763483, 0.9017203154504035, -0.9474931851459008, -0.03334319523220503, -0.14354857449259437, -0.258603947854119, -0.6211074642796408, 0.9328743112042068, 0.5983190378042788, 0.860564215444357, -0.5329857672153024, 0.403783074281117, 0.538582127861995, 0.1061505899839121, -0.9093445419255864, 0.6656150775217203].
    Elapsed time: 1595.9994523525238 s


Example 2: Multi-process Solving
--------------------------------
This example expands on Example 1 and solves the Electron optimization problem with 4 processes.
Note that Python as a language is not great at parallell execution in general, and that spawning processes is quite expensive.
`pygosolnp` spawns processes using a `Multiprocessing`_ `multiprocessing.pool` and shares memory between the processes.
As such, it is easy to cause bugs if you are not familiar with how this works, so please read up on it before using it!
For example, you can not pass lambda functions or class functions to process pools, but global functions should work fine.
If you get different results between single-processing and multi-processing execution your python functions are liekly badly defined for multi-processing.

See full source code on GitHub `/python_examples/example_electron.py`_

.. _`Multiprocessing`: https://docs.python.org/3/library/multiprocessing.html
.. _`/python_examples/example_electron_multiprocessing.py`: https://github.com/KristerSJakobsson/pygosolnp/blob/main/python_examples/example_electron_multiprocessing.py

Below example only changes the previous one in the function call, using 4 processes to run.

.. code-block:: python

        results = pygosolnp.solve(obj_func=obj_func,
                              eq_func=eq_func,
                              eq_values=equality_constraints,
                              par_lower_limit=parameter_lower_bounds,
                              par_upper_limit=parameter_upper_bounds,
                              number_of_restarts=20,
                              number_of_simulations=20000,
                              number_of_processes=4,  # Simulations and processes will be executed in 4 processes
                              seed=443,
                              pysolnp_max_major_iter=100,
                              debug=False)

Running this will yield the output:

::

    Solution 1: 244.1550118432253; Solution 2: 243.9490050190484; Solution 3: 185.78533081425041; Solution 4: 244.07921194485854; Solution 5: 216.19236253370485; Solution 6: 194.1742137471891; Solution 7: 258.6157748268509; Solution 8: 205.72538678938517; Solution 9: 244.0944480181356; Solution 10: 217.4090464122706; Solution 11: 201.58045387715478; Solution 12: 247.70691375326325; Solution 13: 243.92615570955812; Solution 14: 192.3944392661305; Solution 15: 243.93657263760585; Solution 16: 247.17924771908508; Solution 17: 244.06529702108125; Solution 18: 244.29427536763717; Solution 19: 199.69130383979302; Solution 20: 243.99315264179037
    Best solution 243.92615570955812 for parameters [0.8726149386907173, 0.1488320711741995, -0.8215181712229778, 0.8597822831494584, -0.265961670940264, -0.6664127144955102, -0.6029702658967409, 0.2867960203292267, -0.04380531711098636, 0.9519854892760677, -0.39592769694574026, -0.2160514547351913, -0.21416235954836016, 0.4338472533837847, -0.9411378567701716, 0.6418976636970082, 0.014864034847848012, 0.6981416769347426, 0.4413252856284809, -0.5267725521555819, -0.9148568048943023, -0.5831731928212042, 0.47570915153781534, 0.4089885176760918, 0.008471540399374077, -0.36287443863890595, 0.8618964461129363, 0.5476494687199884, -0.3309316231117961, 0.9582851670742292, -0.6505818085537286, 0.2793946112676732, -0.7596998666078645, 0.65142774983249, 0.30572406841664945, -0.1736400992779951, -0.2357569641249718, -0.9762296783338298, 0.8894784482368485, -0.21768032982807542, 0.44966067028074935, 0.359898210796523, 0.3932146838134686, -0.25429503229562933, -0.6621520897149067, 0.0002565729867240561, 0.6081775900274631, -0.8731755460834034, -0.07630776960802095, -0.7462707639808169, 0.32690759610807246, 0.4847543563757037, -0.15870866693945487, -0.38892531575475037, -0.10466177783304143, 0.36421374544164403, -0.7472412325499505, -0.583622807257543, -0.7574487346380878, -0.01614470971763483, 0.9017203154504035, -0.9474931851459008, -0.03334319523220503, -0.14354857449259437, -0.258603947854119, -0.6211074642796408, 0.9328743112042068, 0.5983190378042788, 0.860564215444357, -0.5329857672153024, 0.403783074281117, 0.538582127861995, 0.1061505899839121, -0.9093445419255864, 0.6656150775217203].
    Elapsed time: 596.5835165977478 ms


As expected, we get the same result as in the single-processed approach in example 1!

This is expected as the same random starting values are used as when running single-processed.

Furthermore, the execution time is around one third of the single-processed one. Note that if we reduce the number of restarts to 4 the single-processed execution would be quicker.

Example 3: Truncated Normal Distribution
----------------------------------------
PYGOSOLNP does not depend on any large-scale library (pandas, numpy, scipy etc.) out of box.
This example shows how to overrides the logic for generating starting points by using Scipy and the Truncated Normal distribution.
It is fairly trivial to modify this example to use other `Scipy Distributions`_, for example Beta Distribution sampling etc.

.. _`Scipy Distributions`: https://docs.scipy.org/doc/scipy/reference/stats.html

See full source code on GitHub `/python_examples/example_truncated_normal.py`_

.. _`/python_examples/example_truncated_normal.py`: https://github.com/KristerSJakobsson/pygosolnp/blob/main/python_examples/example_truncated_normal.py

.. code-block:: python

    # The Sampling class is an abstract class that can be inherited and customized as you please
    class TruncatedNormalSampling(pygosolnp.sampling.Sampling):

        def __init__(self,
                     parameter_lower_bounds: List[float],
                     parameter_upper_bounds: List[float],
                     seed: Optional[int]):
            self.__generator = Generator(PCG64(seed))
            self.__parameter_lower_bounds = parameter_lower_bounds
            self.__parameter_upper_bounds = parameter_upper_bounds

        def generate_sample(self, sample_size: int) -> List[float]:
            # This function is abstract, it returns random starting values for one sample
            return truncnorm.rvs(a=self.__parameter_lower_bounds,
                                 b=self.__parameter_upper_bounds,
                                 size=sample_size,
                                 random_state=self.__generator)

    ...
    # See original file for full code
    ...

    # Instantiate sampling object
    sampling = TruncatedNormalSampling(
        parameter_lower_bounds=parameter_lower_bounds,
        parameter_upper_bounds=parameter_upper_bounds,
        seed=99)

    results = pygosolnp.solve(
        obj_func=permutation_function,
        par_lower_limit=parameter_lower_bounds,
        par_upper_limit=parameter_upper_bounds,
        number_of_restarts=6,
        number_of_simulations=2000,
        pysolnp_max_major_iter=25,
        pysolnp_tolerance=1E-9,
        start_guess_sampling=sampling)

Running this will yield:

::

    # Solution 1: 0.0016119745327847497; Solution 2: 0.005968645850086645; Solution 3: 0.006083292803668321; Solution 4: 0.006629107105976147; Solution 5: 0.005305936314073526; Solution 6: 0.006049589559946693
    # Best solution: [1.3008954298086124, 3.181786909056148, 1.3814249752478918, 3.9436695447632877]
    # Objective function value: 0.0016119745327847497
    # Elapsed time: 8.562503099441528 s

With 2000 simulations we got a fairly accurate value in 8.5 seconds.
Lets compare this with Grid Sampling below.

Example 4: Grid Sampling
------------------------
PYGOSOLNP does not depend on any large-scale library (pandas, numpy, scipy etc.) out of box.
This example overrides the logic for generating starting points by using Scikit-optimize Grid Sampling.
It is fairly trivial to modify this example to use other `Scikit-optimize Sampling Methods`_, for example Sobol, Latin hypercube sampling etc.

.. _`Scikit-optimize Sampling Methods`: https://scikit-optimize.github.io/stable/auto_examples/sampler/initial-sampling-method.html

See full source code on GitHub `/python_examples/example_grid_sampling.py`_

.. _`/python_examples/example_grid_sampling.py`: https://github.com/KristerSJakobsson/pygosolnp/blob/main/python_examples/example_grid_sampling.py

.. code-block:: python

    # The Sampling class is an abstract class that can be inherited and customized as you please
    class GridSampling(pygosolnp.sampling.Sampling):

        def __init__(self,
                     parameter_lower_bounds: List[float],
                     parameter_upper_bounds: List[float],
                     seed):
            self.__space = skopt.space.Space(dimensions=zip(parameter_lower_bounds, parameter_upper_bounds))
            self.__seed = seed

        def generate_all_samples(self, number_of_samples: int, sample_size: int) -> List[float]:
            # Overwrite this function to define the behavior when generating starting guesses for all samples
            # By default it calls `generate_sample` number_of_samples time
            grid = skopt.sampler.Grid()
            grid_values = grid.generate(dimensions=self.__space.dimensions,
                                        n_samples=number_of_samples,
                                        random_state=self.__seed)
            return list(chain.from_iterable(grid_values))

        def generate_sample(self, sample_size: int) -> List[float]:
            # This function is abstract
            # Not needed since we are generating a grid for all samples
            pass

    ...
    # See original file for full code
    ...

    # Instantiate sampling object
    sampling = GridSampling(
        parameter_lower_bounds=parameter_lower_bounds,
        parameter_upper_bounds=parameter_upper_bounds,
        seed=92)

    results = pygosolnp.solve(
        obj_func=permutation_function,
        par_lower_limit=parameter_lower_bounds,
        par_upper_limit=parameter_upper_bounds,
        number_of_restarts=6,
        number_of_simulations=2000,
        pysolnp_max_major_iter=25,
        pysolnp_tolerance=1E-9,
        start_guess_sampling=sampling)


Running this will yield the output:

::
    # Solution 1: 0.0006360327708392506; Solution 2: 0.006239163594915304; Solution 3: 0.006140229082904356; Solution 4: 0.006218870214655177; Solution 5: 0.005963823643719209; Solution 6: 0.13065649880545976
    # Best solution: [1.1622677695732497, 1.683172007310748, 3.9509962074974956, 3.159134907203731]
    # Objective function value: 0.0006360327708392506
    # Elapsed time: 22.986207962036133 s

With 2000 simulations Grid Sampling gave a better result than Truncated Normal but it took longer.
