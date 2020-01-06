// @file matrix_alloc.hpp
// @brief  Because I allocate 2-dimensional arrays frequently, I created a 
//         few optional functions that make this more convenient.

#ifndef _MATRIX_ALLOC_H
#define _MATRIX_ALLOC_H

namespace matrix_alloc {

/// @brief  Allocate a 2-dimensional table row-major order
template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// Slightly different version of Alloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
void Alloc2D(size_t M,              //!< size of the array (outer)
             size_t N,              //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// This function is the corresponding way to deallocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry **paX,         //!< pointer to 1-D contiguous-memory array
               Entry ***paaX);      //!< pointer to 2-D multidimensional array



// -------------- IMPLEMENTATION --------------

template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  assert(paX && paaX);

  *paX = new Entry [size[0] * size[1]];

  // Allocate a conventional 2-dimensional
  // pointer-to-a-pointer data structure.
  *paaX = new Entry* [size[1]];
  for(Integer iy=0; iy<size[1]; iy++)
    (*paaX)[iy] = &((*paX)[ iy*size[0] ]);
  // The caller can access the contents of *paX using (*paaX)[i][j] notation.
}

template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (outer)
             size_t ncolumns,       //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  size_t size[2];
  size[0] = ncolumns;
  size[1] = nrows;
  Alloc2D(size, paX, paaX);
}


template<typename Entry>
void Dealloc2D(Entry **paX,          //!< pointer to 1-D contiguous-memory array
               Entry ***paaX)        //!< pointer to 2-D multidimensional array
{
  if (paaX && *paaX) {
    delete [] (*paaX);
    *paaX = nullptr;
  }
  if (paX && *paX) {
    delete [] *paX;
    *paX = nullptr;
  }
}


#endif //#ifndef _MATRIX_ALLOC_H
