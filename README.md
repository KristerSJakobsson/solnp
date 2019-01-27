# C++ SOLNP

WARNING: This source is in development and will not currently run properly.

This is a C++ implementation of the SOLNP algorithm by Yinyu Ye (1989).
The algorithm was originally implemented in Matlab, and have gained some fame through it's R implementation (RSOLNP).
Various implementations of the algorithm exists already, however, this version utilizes the power of DLIB and C++11.

This algorithm solves the general optimization problem on the form:


## Getting Started

Currently, the sources are available as is.
Simply download and use in your project as you please.

### Prerequisites

The CMake file assumes Catch2 is installed globally.

## Running the tests

Simply run the CMake script.

## Built With

* [CMake](https://cmake.org/runningcmake/) - Build tools
* [CLion](https://www.jetbrains.com/clion/) - IDE by JetBrains

## Authors

* **Krister S Jakobsson** - *Implementation*

## License

This project is licensed under the Boost License - see the [license](LICENSE.md) file for details

## Acknowledgments

* **Yinyu Ye** -  Publisher of the original algorithm,
[Original Sources](https://web.stanford.edu/~yyye/matlab/)
* **Rsolnp** - An R implementation of SOLNP, 
[Github repository](https://github.com/cran/Rsolnp)