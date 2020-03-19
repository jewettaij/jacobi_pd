[![Build Status](https://travis-ci.org/jewettaij/jacobi_pd.svg?branch=master)](https://travis-ci.org/jewettaij/jacobi_pd.svg?branch=master)
[![codecov](https://codecov.io/gh/jewettaij/jacobi_pd/branch/master/graph/badge.svg)](https://codecov.io/gh/jewettaij/jacobi_pd)
[![C++11](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://isocpp.org/std/the-standard)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/jacobi_pd)]()
[![License: CC0-1.0](https://licensebuttons.net/p/mark/1.0/88x31.png)](https://creativecommons.org/publicdomain/zero/1.0/)


jacobi_pd
===========

## Description

This repository contains a small C++
[header file](include/jacobi.hpp)
that implements the
[Jacobi eigenvalue algorithm](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm).
It is [free of copyright](https://creativecommons.org/publicdomain/zero/1.0/).

The Jacobi algorithm remains one of the oldest and most popular methods for
diagonalizing dense, square, real, symmetric matrices.

The matrices themselves can be implemented as X\*\* (pointer-to-pointer),
vector\<vector\<X\>\>&, fixed-size arrays,
or any other C or C++ object which supports \[i\]\[j\] indexing.
(Here **X** is any real numeric type.  Complex numbers are not supported.)

*(Memory allocation on the heap is avoided except during instantiation.)*


#### The main feature of this repository is it's [license](LICENSE.md).

As of 2020, no simple *public domain* C++11 code
yet exists for matrix diagonalization.
Other C++ libraries such as Eigen or GSL are typically
much larger and use more restrictive licenses.
*(On several occasions, this has prevented me from including
their code in other open-source projects with incompatible licenses.)*
Some repositories may unwittingly contain code
snippets from other sources, such as
[numerical recipes](http://mingus.as.arizona.edu/~bjw/software/boycottnr.html).
This short repository was written from scratch.
No lines of code were borrowed or adapted from other sources.



*Caveats:* The code in this repository does not run in parallel,
and only works on dense square real symmetric matrices.
However it is reasonably
[short, simple](include/jacobi.hpp), 
[fast](benchmarks/README.md) and
[reliable](.travis.yml).
You can do anything you like with this code.


##  Example usage

```cpp
#include "jacobi.hpp"
using namespace jacobi_public_domain;

// ...

int n = 5;       // Matrix size
double **M;      // A symmetric n x n matrix you want to diagonalize
double *evals;   // Store the eigenvalues here.
double **evects; // Store the eigenvectors here.
// Allocate space for M, evals, and evects, and load contents of M (omitted)...

// Now create an instance of Jacobi ("eigen_calc").  This will allocate space
// for storing intermediate calculations.  Once created, it can be reused
// multiple times without incurring the cost of allocating memory on the heap.

Jacobi<double, double*, double**> eigen_calc(n);

// Note:
// If the matrix you plan to diagonalize (M) is read-only, use this instead:
// Jacobi<double, double*, double**, double const*const*> eigen_calc(n);
// If you prefer using vectors over C-style pointers, this works also:
// Jacobi<double, vector<double>&, vector<vector<double>>&> eigen_calc(n);

// Now, calculate the eigenvalues and eigenvectors of M

eigen_calc.Diagonalize(M, evals, evects);
```

## Benchmarks

[![benchmarks](benchmarks/benchmarks.png)](benchmarks/README.md)

[***(details here...)***](benchmarks/README.md)


## Installation

Copy the file(s) in the [include](include) subdirectory,
to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).
No linking is necessary.
This is a header-only library.


## Development Status: *stable*

**jacobi_pd** has been
[tested](.travis.yml)
for accuracy and memory safety
over a wide range of array types, matrix sizes,
eigenvalue magnitudes and degeneracies.


## Requirements

A C++11 compatible compiler.


## License

*jacobi_pd* is available under the terms of the [Creative-Commons-Zero license](LICENSE.md).

*Please send me corrections or suggestions.*

