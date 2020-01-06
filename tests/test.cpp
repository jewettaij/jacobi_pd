#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <random>
#include "matrix_alloc.hpp"
#include "jacobi.hpp"


using std::cout;
using std::endl;
using std::setprecision;
using namespace matrix_alloc;
using namespace jacobi_public_domain;


template<typename T>
using vector = std::vector<T>;

template<typename T>
using complex = std::complex<T>;

/// @brief
/// Generate a random orthogonal n x n matrix
template<typename Scalar, typename Matrix>
void GenRandOrth(Matrix R,
                 int n,
                 std::default_random_engine &generator=nullptr)
{
  std::normal_distribution<double> gaussian_distribution(0,1);

  for (int i = 0; i < n; i++) {
    // Generate a vector, "v", in a random direction subject to the constraint
    // that it is orthogonal to the first i-1 rows-vectors of the R matrix.
    Scalar rsq = 0.0;
    while (rsq == 0.0) {
      // Generate a vector in a random direction
      // (This works because we are using a normal (Gaussian) distribution)
      for (int j = 0; j < n; j++)
        v[j] = gaussian_distribution(generator);

      //Now subtract from v, the projection of v onto the first i-1 rows of R.
      //This will produce a vector which is orthogonal to these i-1 row-vectors.
      for (int k = 0; k < i; k++) {
        Scalar v_dot_Ri = 0.0
          for (int j = 0; j < n; j++)
            v_dot_Ri += v[j] * R[k][j]; // = <v , R[i]>
        for (int j = 0; j < n; j++)
          v[j] -= v_dot_Ri * R[i][j]; // = v - <V,R[i]> R[i]
      }
      // check if it is linearly independent of the other vectors and non-zero
      rsq = 0.0;
      for (int j = 0; j < n; j++)
        rsq += v[j]*v[j];
    }
    // Now normalize the vector
    r_inv = 1.0 / std::sqrt(rsq);
    for (int j = 0; j < n; j++)
      v[j] *= r_inv;
    // Now copy this vector to the i'th row of R
    for (int j = 0; j < n; j++)
      R[i][j] = v[j];
  } //for (int i = 0; i < n; i++)
} //void GenRandOrth()




/// @brief
/// Multiply two matrices A and B, store the result in C. (C = AB).

template<typename Scalar>
void mmult(Scalar const *const *A, //<! input array
           Scalar const *const *B, //<! input array
           Scalar** C,  //<! store result here
           size_t m,    //<! number of rows of A
           size_t n=0,  //<! optional: number of columns of B (=m by default)
           size_t K=0   //<! optional: number of columns of A = num rows of B
           )
{
  assert(A != C);
  if (n == 0) n = m; // if not specified, assume matrix is square
  if (K == 0) K = m; // if not specified, assume matrix is square

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      C[i][j] = 0.0;

  // perform matrix multiplication
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      for (size_t k = 0; k < K; k++)
        C[i][j] += A[i][k] * B[k][j];
}


template <typename Scalar>
static inline void MatProduct3(Scalar const* const *A,
                               Scalar const* const *B,
                               Scalar **dest) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      dest[i][j] = 0.0;
      for (int k=0; k<3; k++)
        dest[i][j] += A[i][k] * B[k][j];
    }
  }
}



template <typename Scalar>
void TestDiagonalization(int n=3
                         int n_matrices=100,
                         int n_tests_per_matrix=50)
{
  cout << endl << "-- Diagonalization test (real symmetric)  --" << endl;

  Jacobi eigen_calc(n);

  Scalar *evals = new Scalar[n];

  // construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<Scalar> random_real01();

  Scalar **M, **R, **Rt, **D, **tmp;
  Scalar *_M, *_R, *_Rt, *_D, *_tmp;
  Alloc2D(n, n, &_M, &M);
  Alloc2D(n, n, &_R, &R);
  Alloc2D(n, n, &_Rt, &Rt);
  Alloc2D(n, n, &_D, &D);
  Alloc2D(n, n, &_tmp, &tmp);

  for(int imat = 0; imat < n_matrices; imat++) {

    // Randomly generate the eigenvalues
    // D is the diagonal matrix whose diagonal elements are the eigenvalues
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        D[i][j] = 0.0;
    for (int i = 0; i < n; i++) {
      D[i][i] = random_real01();
      if (D[i][i] != 0.0)
        //D[i][i] = logarithm of a random number from (0,1]
        D[i][i] = -std::log(D[i][i]);
      if (random_real01() < 0.5)
        D[i][i] = -D[i][i];
    }

    // Randomly generate the R and Rt matrices:
    GenRandOrth(R, n, generator);
    // Rt is the transpose of R
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        Rt[i][j] = R[j][i];

    // Construct the test matrix M = Rt * D * R
    mmult(Rt, D, tmp, n);
    mmult(tmp, R, M, n);

    for (int i_test = 0; i_test < n_tests_per_matrix; i++) {
      //cout << "Eigenvectors (columns) which are known in advance, Rij = \n";
      //for (int i = 0; i < n; i++) {
      //  for (int j = 0; j < n; j++) {
      //    cout << R[i][j] << " ";
      //  }
      //  cout << "\n";
      //}
      //cout << "  (The eigenvector should match one of these columns.)\n";
      //cout << "\n";

      // Now, calculate the eigenvalues and eigenvectors
      eigen_calc.Jacobi(M, eval, Rt);
      //cout << " --- Results using PEigenCalculator ---\n";
      //cout << "Eigen value: " << setprecision(16) << eval << endl;
      //cout << "Eigen vector:";
      //for(int i = 0;i < n;i++) {
      //  cout << evec[i] << " ";
      //}
      //cout << endl;

    } //for (int i_test = 0; i_test < n_tests_per_matrix; i++)

  } //for(int imat = 0; imat < n_matrices; imat++) {


  
  Dealloc2D(&_M, &M);
  Dealloc2D(&_R, &R);
  Dealloc2D(&_Rt, &Rt);
  Dealloc2D(&_D, &D);
  Dealloc2D(&_tmp, &tmp);
  delete [] evals;

  return EXIT_SUCCESS;
}
