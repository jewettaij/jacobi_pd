/// @file    matrix_alloc.hpp
/// @brief   Because I allocate 2-dimensional arrays frequently, I created a 
///          few optional functions that make this more convenient.
/// @author  Andrew Jewett
/// @license CC0-1.0

#include<cassert>

#ifndef _MATRIX_ALLOC_H
#define _MATRIX_ALLOC_H

namespace matrix_alloc {

/// @brief Allocate a 2-dimensional table row-major order.
template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (number of rows)
             size_t ncols,          //!< size of the array (number of columns)
             Entry ***paaX          //!< pointer to a 2-D C-style array
             );

/// @brief  This function is the corresponding way to deallocate
///         arrays that were created using Alloc2D().
template<typename Entry>
void Dealloc2D(Entry ***paaX        //!< pointer to 2-D multidimensional array
               );


// --- IMPLEMENTATION ---

template<typename Entry>
void Alloc2D(size_t nrows,          // size of the array (number of rows)
             size_t ncols,          // size of the array (number of columns)
             Entry ***paaX)         // pointer to a 2-D C-style array
{
  assert(paaX);
  // Allocate a conventional 2-dimensional
  // pointer-to-a-pointer data structure.
  *paaX = new Entry* [nrows];
  (*paaX)[0] = new Entry [nrows * ncols];
  for(size_t iy=0; iy<nrows; iy++)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
  // The caller can access the contents using (*paaX)[i][j] notation.
}

template<typename Entry>
void Dealloc2D(Entry ***paaX)       // pointer to 2-D multidimensional array
{
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}

} // namespace matrix_alloc

#endif //#ifndef _MATRIX_ALLOC_H
