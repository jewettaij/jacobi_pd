// @file matrix_alloc.hpp
// @brief  Because I allocate 2-dimensional arrays frequently, I created a 
//         few optional functions that make this more convenient.

#include<cassert>

#ifndef _MATRIX_ALLOC_H
#define _MATRIX_ALLOC_H

namespace matrix_alloc {

/// @brief  Allocate a 2-dimensional table row-major order
template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// Slightly different version of Alloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
void Alloc2D(size_t M,              //!< size of the array (number of rows)
             size_t N,              //!< size of the array (number of columns)
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// This function is the corresponding way to deallocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry ***paaX);      //!< pointer to 2-D multidimensional array



// -------------- IMPLEMENTATION --------------

template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (number of rows)
             size_t ncolumns,       //!< size of the array (number of columns)
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  assert(paaX);

  Entry *aX = new Entry [size[0] * size[1]];

  // Allocate a conventional 2-dimensional
  // pointer-to-a-pointer data structure.
  *paaX = new Entry* [size[1]];
  for(Integer iy=0; iy<size[1]; iy++)
    (*paaX)[iy] = &(aX[ iy*size[0] ]);
  // The caller can access the contents of *paX using (*paaX)[i][j] notation.
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
