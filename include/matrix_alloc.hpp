// @file matrix_alloc.hpp
// @brief  Because I allocate 2-dimensional arrays frequently, I created a 
//         few optional functions that make this more convenient.

#include<cassert>

#ifndef _MATRIX_ALLOC_H
#define _MATRIX_ALLOC_H

namespace matrix_alloc {

/// @brief
/// Allocate a 2-dimensional table row-major order.
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
void Alloc2D(int M,                 //!< size of the array (number of rows)
             int N,                 //!< size of the array (number of columns)
             Entry ***paaX);        //!< pointer to 2-D multidimensional array


/// @brief Allocate a 2-dimensional table.  After allocation, the entries in
///        the table (*paaX)[i][j] exist for all 0<=i<size[0] and 0<=j<size[1].
template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// This function is the corresponding way to deallocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry ***paaX);      //!< pointer to 2-D multidimensional array


// --- IMPLEMENTATION ---

template<typename Entry>
void Alloc2D(int nrows,             //!< size of the array (outer)
             int ncolumns,          //!< size of the array (inner)
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  int size[2];
  size[0] = nrows;
  size[1] = ncolumns;
  Alloc2D<Entry, int>(size, paaX);
}

template<typename Entry>
void Dealloc2D(Entry ***paaX)        //!< pointer to 2-D multidimensional array
{
  if (paaX && *paaX) {
    Entry *aX = &((*paaX)[0][0]);
    delete [] (*paaX);
    delete [] aX;
    *paaX = nullptr;
  }
}

} // namespace matrix_alloc

#endif //#ifndef _MATRIX_ALLOC_H
