[![Build Status](https://travis-ci.org/jewettaij/jacobi.svg?branch=master)](https://travis-ci.org/jewettaij/jacobi.svg?branch=master)
[![GitHub](https://img.shields.io/github/license/jewettaij/jacobi)](./LICENSE.md)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/jacobi)]()


jacobi
===========

## Description

This repository contains
[public domain](https://creativecommons.org/publicdomain/zero/1.0/)
header-only template C++ source code for the
[Jacobi eigenvalue algorithm](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm).
This remains one of the oldest and most popular algorithms for
diagonalizing dense, square, real, symmetric matrices.

The matrices themselves can be implemented as \*\*X (pointer-to-pointer),
vector\<vector\<X\>\>, fixed-size arrays,
or any other C or C++ object which supports double-indexing.
(Where **X** is any real numeric type.  Complex numbers are not supported.)

*(Memory allocation on the heap is avoided except during initialization.)*


#### The main feature of this repository is it's [license](LICENSE.md).

I have recently been having a hard time finding a short, simple eigenvector calculator in C++ with an explicitly-stated permissive open-source license.
Amazingly, in early 2020, no *public-domain*
C++ code yet exists for matrix diagonalization.
Other C++ libraries such as Eigen or GSL
use somewhat more restrictive licenses.
*(On several occasions, those licenses have prevented me from borrowing code
from those libraries to contribute to other open-source projects.)*
Out of that frustration, I wrote this repository.

*Caveats:* The code in this repository has not been optimized,
does not run in parallel,
and it only works on dense square real symmetric matrices.
However you can freely use this humble code anywhere you like.
Use at your own risk.


## Installation

This is a header-only library.

Copy the files in the [include](include) subdirectory,
to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).

# Development Status: *alpha*

As of 2020-1-23, basic functionality appears to be working.
More testing is needed including tests for memory leaks and profiling.
The API might change slightly, but existing code built using
it should still work.
(Later I might add "const" to the first argument of Jacobi::Diagonalize().)


## License

*jacobi* is available under the terms of the [Creative-Commons-Zero license](LICENSE.md).

*Please send me corrections or suggestions.*

