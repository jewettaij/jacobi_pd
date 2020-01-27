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
(Here **X** is any real numeric type.  Complex numbers are not supported.)

*(Memory allocation on the heap is avoided except during initialization.)*


#### The main feature of this repository is it's [license](LICENSE.md).

I have recently been having a hard time finding a short, simple eigenvector
calculator in C++ with an explicitly-stated permissive open-source license.
Amazingly, in early 2020, no simple *public-domain*
C++ code yet exists for matrix diagonalization.
Other C++ libraries such as Eigen or GSL use somewhat more restrictive licenses.
*(On several occasions, those licenses have prevented me from borrowing code
from these libraries to contribute to other open-source projects.)*  Some
repositories may unwittingly contain code snippets from other sources, such as
[numerical recipes](http://mingus.as.arizona.edu/~bjw/software/boycottnr.html).
This repository was written from scratch.  No lines of code were borrowed
or adapted from other sources.



*Caveats:* The code in this repository has not been optimized,
does not run in parallel,
and it only works on dense square real symmetric matrices.
However you can freely use this code anywhere you like.

##  Example usage

```cpp
#include "jacobi.hpp"
using namespace jacobi_public_domain;

// ...

double **M;      // A symmetric matrix you want to diagonalize
double *evals;   // Store the eigenvalues here.
double **evects; // Store the eigenvectors here.
// Allocate space for M, evals, and evects, and load contents of M (omitted)...

// Now create an instance of Jacobi ("eigen_calc").  This will allocate space
// for storing intermediate calculations.  Once created, it can be reused
// multiple times without incurring the cost of allocating memory on the heap.
int n = 5;
Jacobi<double, double*, double**> eigen_calc(n);

// Now, calculate the eigenvalues and eigenvectors of M
eigen_calc.Diagonalize(M, evals, evects);
```
*(A complete working example can be found [here](tests/test.cpp).)*

## Installation

Copy the file(s) in the [include](include) subdirectory,
to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).
No linking is necessary.
This is a header-only library.

# Development Status: *alpha*

As of 2020-1-23, basic functionality appears to be working.
More testing is needed including tests for memory leaks and profiling.
The API might change slightly, but existing code built using
it should still work.
(Later I might add "const" to the first argument of Jacobi::Diagonalize().)


## License

*jacobi* is available under the terms of the [Creative-Commons-Zero license](LICENSE.md).

*Please send me corrections or suggestions.*

