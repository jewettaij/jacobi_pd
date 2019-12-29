[![Build Status](https://travis-ci.org/jewettaij/jacobi.svg?branch=master)](https://travis-ci.org/jewettaij/jacobi.svg?branch=master)
[![GitHub](https://img.shields.io/github/license/jewettaij/jacobi)](./LICENSE.md)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/jacobi)]()


jacobi
===========

This repository contains [***public-domain***](LICENSE.md)
header-only C++ source code for the
[Jacobi eigenvalue algorithm](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm).
This remains one of the oldest and most popular algorithms for
diagonalizing small, dense, square, real, symmetric matrices.
The matrices themselves can be implemented as \*\*X (pointer-to-pointer),
vector\<vector\<X\> \>, fixed-size arrays,
or any other C++ object which supports double-indexing.
(**X** is any real numeric type.  Complex numbers are not supported.)

#### The main feature of this repository is it's license.

Amazingly, as of 2020-1-01, no public-domain C++ header code exists for matrix diagonalization.  Other C++ libraries such as Eigen or GSL use more restrictive licenses.  (In my case, that prevent me from using those libraries to
contribute to some github projects.  So I wrote my own library.)
This code has not been optimized, and does not run in parallel,
and it only works on dense square real symmetric matrices.
Use at your own risk.

## Development Status: *Planning (pre-alpha)*

As of 2019-12-28, this code is not complete and does not compile.

## Installation

This is a header-only library.

Copy the files in the [include/jacobi.hpp](include) subdirectory,
to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).

## License

*jacobi* is available under the terms of the [Creative-Commons-Zero](LICENSE.md).

*Please send me corrections or suggestions.*

