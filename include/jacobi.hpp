///   @file  jacobi.hpp
///   @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
///          using the Jacobi eigenvalue algorithm.
#ifndef _JACOBI_HPP
#define _JACOBI_HPP

#include <cmath>
#include <cassert>


namespace jacobi {

// Because I allocate 2-dimensional arrays frequently, I created a 
// few functions that make this more convenient.

/// @brief  Allocate a 2-dimensional table row-major order
template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// Slightly different version of Alloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
void Alloc2D(int M,                 //!< size of the array (outer)
             int N,                 //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// This function is the corresponding way to dellocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry **paX,          //!< pointer to 1-D contiguous-memory array
               Entry ***paaX);       //!< pointer to 2-D multidimensional array


template<typename Scalar>
static inline Scalar SQR(Scalar x) {return x*x;}


/// @class Jacobi
/// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
///        using the Jacobi eigenvalue algorithm.
///        The algorithm implemented here follows the strategy explained here:
///        https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
///        https://web.stanford.edu/class/cme335/lecture7.pdf
/// @note  The "Vector", "Matrix" and "ConstMat" type arguments can be any 
///        C or C++ object that support indexing, including pointers or vectors.

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
class Jacobi
{
  int n;            // the size of the matrix
  // The next 3 data members store the rotation, translation and scale
  // after optimal superposition
  Scalar **M;              //!< store copy of matrix here (workspace)
  Scalar *_M;              //!< contents of M (contiguously allocated)
  int *max_ind_per_row;    //!< for each row, store the index of the max entry
  bool   *changed_row;     //!< was this row changed during previous iteration?
  // Precomputed cosine, sin, and tangent of the most recent rotation angle:
  Scalar c;                //!< = cos(θ)
  Scalar s;                //!< = sin(θ)
  Scalar t;                //!< = tan(θ),  (note |t|<=1)


public:
  void SetSize(int matrix_size) {
    n = matrix_size;
    Dealloc();
    Alloc(n);
  }

  int GetSize() {
    return n;
  }

  Jacobi(int matrix_size=0) {
    Init();
    SetSize(matrix_size);
  }

  /// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
  ///        using the Jacobi eigenvalue algorithm:
  ///        https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
  void
  Diagonalize(ConstMat M,     //!< the matrix you wish to diagonalize (size n)
              Vector eval,    //!< store the eigenvalues here
              Matrix evec);   //!< store the eigenvectors here (in rows)

private:

  /// @brief Calculate the components of a rotation matrix which performs a
  ///        i,j plane by an angle (θ) that (when multiplied on both sides)
  ///        will zero the ij'th element of M, so that afterwards M[i][j] = 0
  ///        (This will also zero M[j][i] since M is assumed to be symmetric.)
  ///        The results will be stored in the c, s, and t data members
  ///        (c,s,t store cos(θ), sin(θ), and tan(θ), respectively.)
  void CalcRot(Matrix M,    //!< matrix
               int i,       //!< row index
               int j)       //!< column index

  // memory management:
  void Alloc(int N);
  void Init();
  void Dealloc();
  // memory management: copy constructor, swap, and assignment operator
  Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMat>& source);
  void swap(Jacobi<Scalar, Vector, Matrix, ConstMat> &other);
  Jacobi<Scalar, Vector, Matrix, ConstMat>& operator = (Jacobi<Scalar, Vector, Matrix, ConstMat> source);

}; // class Jacobi





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
void Alloc2D(int nrows,          //!< size of the array (outer)
             int ncolumns,       //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  int size[2];
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


/// @brief Calculate the components of a rotation matrix which performs a
///        i,j plane by an angle (θ) that (when multiplied on both sides)
///        will zero the ij'th element of M, so that afterwards M[i][j] = 0
///        (This will also zero M[j][i] since M is assumed to be symmetric.)

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
CalcRot(Matrix M,    //!< matrix
        int i,       //!< row index
        int j)       //!< column index
{
  t = 1.0; // = tan(θ)
  Scalar M_jj_ii = (M[j][j] - M[i][i]);
  if (M_jj_ii != 0.0) {
    // kappa = (M[j][j] - M[i][i]) / (2*M[i][j])
    Scalar kappa = M_jj_ii;
    t = 0.0;
    M_ij = M[i][j];
    if (M_ij != 0.0) {
      kappa /= (2.0*M_ij);
      // t satisfies: t^2 + 2*t*kappa - 1 = 0
      // -->  t = -kappa +/- sqrt(1+kappa^2)
      // (choose the root which has the smaller absolute value)
      t = std::sqrt(1 + kappa*kappa);
      if (kappa < 0.0)
        t = -t;
      t -= kappa;
    }
  }
  assert(std::abs(t) <= 1.0);
  c = 1.0 / std::sqrt(1 + t*t);
  s = c*t;
}


/// @brief  Perform a similarity transform by multiplying matrix M on both
///         sides by a rotation matrix (applying its inverse to the other side)
///         This matrix performs a rotation in the i,j plane by angle θ.
/// @details This function assumes that i<j and that cos(θ), sin(θ), and tan(θ)
///         have already been computed (by invoking CalcRot()).
///         To save time, since the matrix is symmetric, the elements
///         below the diagonal (ie. M[u][v] where u>v) are not computed.
/// @code
///   M' = R^T * M * R
/// where R the rotation in the i,j plane and ^T denotes the transpose.
///                 i         j
///       _                             _
///      |  1                            | 
///      |    .                          |
///      |      .                        |
///      |        .                      |
///      |          c   ...   s          |
///      |          .  .      .          |
/// R  = |          .    1    .          |
///      |          .      .  .          |
///      |          -s  ...   c          |
///      |                      .        |
///      |                        .      |
///      |                          .    |
///      |_                           1 _|
///
/// Let M' denote the matrix M after multiplication by R^T and R.
/// The components of M' are:
///   M'_uv =  Σ_w  Σ_z   R_wu * M_wz * R_zv
/// @endcode

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
ApplyRot(Matrix M,  //!< matrix
         int i,     //!< row index
         int j)     //!< column index
{
  // Recall that:
  // c = cos(θ)
  // s = sin(θ)
  // t = tan(θ) (which should be <= 1.0)

  // Compute the diagonal elements of M which have changed:
  assert(std::abs(t) <= 1.0);
  M[i][i] -= t * M[i][j];
  M[j][j] += t * M[i][j];
  // Note: This is algebraically equivalent to:
  // M[i][i] = c*c*M[i][i] + s*s*M[j][j] - 2*s*c*M[i][j]
  // M[j][j] = s*s*M[i][i] + c*c*M[j][j] + 2*s*c*M[i][j]

  // Compute the off-diagonal elements of M which have changed:
  //   This code is slow, but easier to undestand.  Commenting out:
  //for (int w = 0; w < n; w++) {
  //  if (w != i)
  //    M[w][i] = c*M[w][i] - s*M[w][j];
  //  if (w != j)
  //    M[w][j] = s*M[w][i] + c*M[w][j];
  //  M[i][w] = M[w][i];
  //  M[j][w] = M[w][j];
  //}
  // Instead, to save time, we avoid the elements below the diagonal

  assert(i < j);
  // compute M[w][i] and M[i][w] for all w!=i
  for (int w=0; w < i; w++)
    M[w][i] = c*M[w][i] - s*M[w][j];    // 0 <= w < i <  j < n
  for (int w=i+1; w < j; w++)
    M[i][w] = c*M[i][w] - s*M[w][j];    // 0 <= i < w <  j < n
  for (int w=j; w < n; w++)
    M[i][w] = c*M[i][w] - s*M[j][w];    // 0 <= i < j <= w < n

  // compute M[w][j] and M[w][j] for all w!=j
  for (int w=0; w < i; w++)
    M[w][j] = s*M[w][i] + c*M[w][j];    // 0 <= w <  i <  j < n
  for (int w=i; w < j; w++)
    M[w][j] = s*M[i][w] + c*M[w][j];    // 0 <= i <= w <  j < n
  for (int w=j+1; w < n; w++)
    M[j][w] = s*M[i][w] + c*M[j][w];    // 0 <= i <  j <  w < n
}


/// @brief  Multiply matrix M on the RIGHT side only by a rotation matrix.
///         This matrix performs a rotation in the i,j plane by angle θ  (where
///         the arguments "s" and "c" refer to cos(θ) and sin(θ), respectively).
/// @code
///
/// E'_uv = Σ_w  E_uw * R_wv
///
/// where:
///                 i         j
///       _                             _
///      |  1                            | 
///      |    .                          |
///      |      .                        |
///      |        .                      |
///      |          c   ...   s          |
///      |          .  .      .          |
/// R  = |          .    1    .          |
///      |          .      .  .          |
///      |          -s  ...   c          |
///      |                      .        |
///      |                        .      |
///      |                          .    |
///      |_                           1 _|
/// @endcode

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
ApplyRotRight(Matrix E,  //!< matrix
              int i,     //!< row index
              int j)     //!< column index
{
  // Recall that c = cos(θ) and s = sin(θ)
  for (int u = 0; u < n; u++) {
    E[u][i] = c*E[u][i] - s*E[u][j];
    E[u][j] = s*E[u][i] + c*E[u][j];
  }
}



template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
MaxEntryPerRow(ConstMat M, int i) const {
  int j_max = i+1;
  for(int j = i+2; j < n)
    if (std::abs(M[i][j]) > std::abs(M[i][j_max]))
      j_max = j;
  assert(j_max > i);
  return j_max;
}



template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
MaxEntry(ConstMat M, int& i_max, int& j_max) const {
  // find the maximum entry in the matrix M in O(n) time
  Scalar max_entry = 0.0;
  i_max = 0;
  j_max = 0;
  for (int i=0; i < n; i++) {
    j = max_entry_per_row[i];
    if (std::abs(M[i][j]) > max_entry) {
      max_entry = std::abs(M[i][j]);
      i_max = i;
      j_max = j;
    }
  }
  #ifndef NDEBUG
  // -- remove the next 4 lines before publishing --
  // make sure that the maximum element really is stored at i_max, j_max
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      assert(std::abs(M[i][j]) <= max_entry);
  // --
  #endif
}


template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
UpdateMaxEntry(ConstMat M, int i, int j) {
  // Update "max_entry_per_row[]" after a Givens rotation at location i,j.
  // Note that a Givens rotation at location i,j will modify all of the matrix
  // elements containing at least one index which is either i or j
  // such as: M[w][i], M[i][w], M[w][j], M[j][w].
  // Check and see whether these modified matrix elements exceed the 
  // corresponding values in max_entry_per_row[] array for that row.
  // If so, then update max_entry_per_row for that row.
  // This is somewhat complicated by the fact that we must only consider
  // matrix elements in the upper-right triangle
  // (ie. matrix elements whose second index is >= the first index).
  // The modified elements we must consider are marked with an "X" below:
  //                 i         j
  //       _                             _
  //      |  .       X         X          | 
  //      |    .     X         X          |
  //      |      .   X         X          |
  //      |        . X         X          |
  //      |          X X X X X X X X X X  |
  //      |            .       X          |
  //      |              .     X          |
  // M  = |                .   X          |
  //      |                  . X          |
  //      |                    X X X X X  |
  //      |                      .        |
  //      |                        .      |
  //      |                          .    |
  //      |_                           . _|

  Scalar tmp;
  for (int w = 0; w < i; w++) {
    //M[w][i] was modified.  See if it exceeds the max_entry on row w
    tmp = std::abs(M[w][i]);
    if (tmp > max_entry_per_row[w])
      max_entry_per_row[w] = i;
    //M[w][j] was modified.  See if it exceeds the max_entry on row w
    tmp = std::abs(M[w][j]);
    if (tmp > max_entry_per_row[w])
      max_entry_per_row[w] = j;
  }
  for (int w = i+1; w < j; w++) {
    //M[w][j] was modified.  See if it exceeds the max_entry on row w
    tmp = std::abs(M[w][j]);
    if (tmp > max_entry_per_row[w])
      max_entry_per_row[w] = j;
  }
  for (int w = i; w < n; w++) {
    //M[i][w] was modified.  See if it exceeds the max_entry on row i
    tmp = std::abs(M[i][w]);
    if (tmp > max_entry_per_row[i])
      max_entry_per_row[i] = w;
    if (w >= j) {
    //M[j][w] was modified.  See if it exceeds the max_entry on row j
    tmp = std::abs(M[j][w]);
    if (tmp > max_entry_per_row[j])
      max_entry_per_row[j] = w;
    }
  }
}



template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
Diagonalize(ConstMat M,     //!< the matrix you wish to diagonalize (size n)
            Vector eval,    //!< store the eigenvalues here
            Matrix evec)    //!< store the eigenvectors here (in rows)
{

  while (NOT_CONVERGED) {
    int i,j;
    for (
  }
}






template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
Init() {
  n = 0;
  M = nullptr;
  max_ind_per_row = nullptr;
  changed_row = nullptr;
}

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
Alloc(int n) {
  this->n = n;
  Alloc2D(n, n, &_M, &M);
  max_ind_per_row = new int[n];
  changed_row = new bool[n];
  assert(_M && M && max_ind_per_row && changed_row);
}

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
Dealloc() {
  if (M)
    Dealloc2D(&_M, &M);
  if (max_ind_per_row)
    delete [] max_ind_per_row;
  if (changed_row)
    delete [] changed_row;
  Init();
}

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
Jacobi<Scalar, Vector, Matrix, ConstMat>::
Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMat>& source)
{
  Init();
  Alloc(source.n);
  assert(n == source.n);
  std::copy(source._M
            source._M + n*n,
            _M);
  std::copy(source.max_ind_per_row,
            source.max_ind_per_row + n,
            max_ind_per_row);
  std::copy(source.changed_row,
            source.changed_row + n,
            changed_row);
}

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
void Jacobi<Scalar, Vector, Matrix, ConstMat>::
swap(Jacobi<Scalar, Vector, Matrix, ConstMat> &other) {
  std::swap(_M, other._M);
  std::swap(max_ind_per_atom, other.max_ind_per_atom);
  std::swap(changed_row, other.changed_row);
  std::swap(n, other.n);
}

template<typename Scalar, typename Vector, typename Matrix, typename ConstMat>
Jacobi<Scalar, Vector, Matrix, ConstMat>&
Jacobi<Scalar, Vector, Matrix, ConstMat>::
operator = (Jacobi<Scalar, Vector, Matrix, ConstMat>y source) {
  this->swap(source);
  return *this;
}



} // namespace jacobi


#endif //#ifndef _JACOBI_HPP

