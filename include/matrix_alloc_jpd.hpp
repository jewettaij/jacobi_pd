/// @file    matrix_alloc_jpd.hpp
/// @brief   Simple functions for allocating 2-dimensional C-style ** arrays in
///          contiguous memory. Perhaps it is useful to put this short code in a
///          separate file to make it available independently of "jacobi_pd.hpp"
/// @author  Andrew Jewett
/// @license CC0-1.0

#ifndef MATRIX_ALLOC_JPD_HPP_
#define MATRIX_ALLOC_JPD_HPP_

#include<cassert>

namespace matrix_alloc_jpd {

/// @brief  Allocate a 2-dimensional array.  (Uses row-major order.)
/// @param  nrows  size of the array (number of rows)
/// @param  nrows  size of the array (number of columns)
/// @param  paaX   pointer to a 2D C-style array
template<typename Entry>
void Alloc2D(std::size_t nrows, std::size_t ncols, Entry ***paaX);

/// @brief  Deallocate arrays that were created using Alloc2D().
/// @param  paaX   pointer to a 2D C-style array
template<typename Entry>
void Dealloc2D(Entry ***paaX);



// ---- IMPLEMENTATION ----


template<typename Entry>
void Alloc2D(std::size_t nrows, std::size_t ncols, Entry ***paaX) {
  assert(paaX);
  *paaX = new Entry* [nrows];  // conventional 2D C array (pointer-to-pointer)
  (*paaX)[0] = new Entry[nrows * ncols];  // 1D C array (contiguous memory)
  for (std::size_t iy = 0; iy < nrows; ++iy)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
  // The caller can access the contents using (*paaX)[i][j]
}

template<typename Entry>
void Dealloc2D(Entry ***paaX) {
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}

}  // namespace matrix_alloc_jpd

#endif  // MATRIX_ALLOC_JPD_HPP_
