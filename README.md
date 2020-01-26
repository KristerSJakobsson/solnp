# Python/C++ SOLNP

This is a C++ implementation of the SOLNP algorithm by Yinyu Ye (1989) with Python Wrappers.
The algorithm was originally implemented in Matlab, and have gained some fame through it's R implementation (RSOLNP).
Various implementations of the algorithm exists already, however, this version utilizes the power of DLIB and C++11.

This algorithm solves the general nonlinear optimization problem on the form:
```
    minimize f(x)
      subject to
       g(x) = 0
   l_h <= h(x) <= u_x
   l_x <=  x   <= u_X
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

The sources for all prerequisites are downloaded by the CMake script when run.
These are:
- pybind11 - Bindings for C++ to Python
- dlib - A c++ mathematical library
- catch2 (for tests) - A testing library

## Running the tests

Run the tests `solnp_test` using CMake.

## Built With

* [CMake](https://cmake.org/runningcmake/) - Build tools
* [CLion](https://www.jetbrains.com/clion/) - IDE by JetBrains

## Authors

* **Krister S Jakobsson** - *Implementation* - krister.s.jakobsson@gmail.com

## License

This project is licensed under the Boost License - see the [license](LICENSE.md) file for details

## Acknowledgments

* **Yinyu Ye** -  Publisher of the original algorithm,
[Original Sources](https://web.stanford.edu/~yyye/matlab/)
* **Rsolnp** - An R implementation of SOLNP, 
[Github repository](https://github.com/cran/Rsolnp)