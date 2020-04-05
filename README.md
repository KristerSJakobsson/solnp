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

## CI and building Wheels

This project uses CI to automatically build wheels for a wide range of distributions.
Notably currently only builds for CPython are available on PyPi, but one can also manually installing the package from source as explained above.

Appveyor CI:
  - Windows with Visual Studio - PARTIALLY IMPLEMENTED
  - Mac OS with CMake - NOT IMPLEMENTED
  
Travis CI:
  - `manylinux2014` Docker with GCC

## Built With

Libraries:
* [dlib](http://dlib.net/) - C++ math library
* [pybind11](https://github.com/pybind/pybind11) - Bindings for building Python Wheels with c++11
* [manylinux](https://github.com/pypa/manylinux) - Docker images for building Linux wheels

Tools:
* [CMake](https://cmake.org/runningcmake/) - Build tools
* [CLion](https://www.jetbrains.com/clion/) - IDE by JetBrains
* [Travis CI](https://travis-ci.org/) - Travis CI for building Manylinux Wheels
* [Appveyor CI](https://www.appveyor.com/) - Appveyor CI for building Windows and Mac OS Wheels

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
