[![codecov](https://codecov.io/gh/KristerSJakobsson/solnp/branch/master/graph/badge.svg)](https://codecov.io/gh/KristerSJakobsson/solnp)
[![Documentation Status](https://readthedocs.org/projects/solnp/badge/?version=latest)](https://solnp.readthedocs.io/en/latest/?badge=latest)

See full documentation on [http://solnp.readthedocs.io](https://solnp.readthedocs.io/en/latest/).

# Python/C++ SOLNP

This is a C++ implementation of the SOLNP algorithm by Yinyu Ye (1989) with Python Wrappers.
The algorithm was originally implemented in Matlab, and have gained some fame through it's R implementation (RSOLNP).
Various implementations of the algorithm exists already, however, this version utilizes the power of DLIB and C++11.

This algorithm solves the general nonlinear optimization problem on the form:
```
    minimize f(x)
      subject to
       g(x) = 0
   l_h <= h(x) <= u_h
   l_x <   x   < u_X
```
where f(x), g(x) and h(x) are smooth functions.

## Getting Started in Python

Simply install the package:
`pip install pysolnp`
<br>
See the `/python_examples` folder for examples.

## Getting Started in C++

The files are header-only and only rely on dlib.
Import `solnp.hpp` and call the solnp function.
<br>
See the `/test` folder for examples.

### Prerequisites

The sources for all prerequisites are linked using github submodules.
To compile the tests, run the CMake script.

Prerequisites for running the C++ SOLNP tests are:
- dlib - A C++ mathematical library
- catch2 (for tests) - A testing library

Additionally, when building the Python wheels you need:
- pybind11 - Bindings for C++ to Python

## Running the tests

Run the tests in the `solnp_test` and export the results to xml by running below in the solnp root folder:
```
$ cmake .
$ make solnp_tests
$ ./solnp_tests -r junit > solnp_tests_result.xml
$ make utils_tests
$ ./utils_tests -r junit > utils_tests_result.xml
```


## CI and building Wheels

This project uses CI to automatically build wheels for a wide range of distributions.
Notably currently only builds for CPython are available on PyPi, but one can also manually installing the package from source as explained above.
Apple Silicon (M1 etc) compiling is not currently available on Github Actions. 

Github Actions:
  - Windows with Visual Studio
  - Mac OS with Clang
  - `manylinux2014` Docker with GCC
  - CodeCov

ReadTheDocs CI:
  - Building and hosting the [docs](https://solnp.readthedocs.io/en/latest/)

## Built With

Libraries:
* [dlib](http://dlib.net/) - C++ math library
* [pybind11](https://github.com/pybind/pybind11) - Bindings for building Python Wheels with C++11
* [manylinux](https://github.com/pypa/manylinux) - Docker images for building Linux wheels
* [cibuildwheels](https://cibuildwheel.readthedocs.io/en/stable) - Library for building Python Wheels through CI

Tools:
* [CMake](https://cmake.org/runningcmake/) - Build tools
* [CLion](https://www.jetbrains.com/clion/) - IDE by JetBrains
* [Github Actions](https://github.com/features/actions) - Building binary Wheels

## Authors

* **Krister S Jakobsson** - *Implementation* - krister.s.jakobsson@gmail.com

## License

This project is licensed under the Boost License - see the [license](LICENSE.md) file for details

## Acknowledgments

* **Yinyu Ye** -  Publisher and mastermind behind the original SOLNP algorithm,
[Original Sources](https://web.stanford.edu/~yyye/matlab/)
* **Alexios Ghalanos and Stefan Theussl** - The people behind RSOLNP,
[Github repository](https://github.com/cran/Rsolnp)
* **Davis King** - The mastermind behind Dlib, check out his blog! [Blog](http://blog.dlib.net/) 
