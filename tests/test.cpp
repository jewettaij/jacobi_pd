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

// are two numbers "similar"?
template<typename T>
inline static bool Similar(T a, T b, T eps=1.0e-06) {
  return std::abs(a - b) < std::abs(eps);
}

// are two vectors (containing n numbers) similar?
template<typename T>
inline static bool Similar(T a, T b, int n, T eps=1.0e-06) {
  for (int i = 0; i < n; i++)
    if (not Similar(a[i], b[i], eps))
      return false;
  return true;
}

//Sort the rows of a matrix "M" by the numbers contained in "keys" (also sorted)
//(This is a simple O(n^2) sorting method, but O(n^2) is a lower bound anyway.)
template<typename Vector, typename Matrix>
void SortRows(Vector keys, Matrix M, int n)
{
  for (int i = 0; i < n; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++)
      if (keys[j] > keys[i_max])
        i_max = j;
    std::swap(keys[i], keys[i_max]);
    for (int k = 0; k < n; k++)
      std::swap(M[i][k], M[i_max][k]);
  }
}


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
void TestJacobi(int n, //<! matrix size
                int n_matrices=100, //<! number of matrices to test
                Scalar eval_magnitude_range=1.0, //<! range of eigevalues
                int n_tests_per_matrix=50) //<! repeat test for benchmarking?
{
  cout << endl << "-- Diagonalization test (real symmetric)  --" << endl;

  // Convert from base 10 to base e and (divide by 2)
  eval_magnitude_range *= std::log(10) * 0.5;

  Jacobi eigen_calc(n);

  Scalar *evals_known = new Scalar[n];
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

    // D is a diagonal matrix whose diagonal elements are the eigenvalues
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        D[i][j] = 0.0;
    for (int i = 0; i < n; i++) {
      // generate some random eigenvalues
      // Use a "log-uniform distribution" (a.k.a. "reciprocal distribution")
      // (This is a way to specify numbers with a precise range of magnitudes.)
      Scalar rand_erange = (random_real01()-0.5); // a number between [-0.5,0.5]
      rand_erange *= eval_magnitude_range; // scale this by eval_magnitude_range
      evals_known[i] = std::exp(rand_erange);
      // also consider both positive and negative eigenvalues:
      if (random_real01() < 0.5)
        evals_known[i] = -evals_known[i];
      D[i][i] = evals_known[i][i];
    }

    // Now randomly generate the R and Rt matrices (ie. the eigenvectors):
    GenRandOrth(R, n, generator);
    // Rt is the transpose of R
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        Rt[i][j] = R[j][i];

    // Construct the test matrix, M, where M = Rt * D * R
    mmult(Rt, D, tmp, n);
    mmult(tmp, R, M, n);

    // Sort the matrix evals and eigenvector rows
    SortRows(evals, R, n);

    if (n_matrices == 1) {
      cout << "Eigenvalues (after sorting):\n"
      for (int i = 0; i < n; i++)
        cout << evals[i] << " ";
      cout << "\n";
      cout << "Eigenvectors (rows) which are known in advance:\n";
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
          cout << R[i][j] << " ";
        cout << "\n";
      }
      cout << "  (The eigenvectors calculated by Jacobi::Diagonalize() should match these.)\n";
      cout << "\n";
    }

    for (int i_test = 0; i_test < n_tests_per_matrix; i++) {

      // Now, calculate the eigenvalues and eigenvectors
      eigen_calc.Jacobi(M, evals, Rt);

      if ((n_matrices == 1) && (i_test == 0)) {
        cout << "Eigenvectors calculated by Jacobi::Diagonalize()\n";
        for (int i = 0; i < n; i++)
          cout << evals[i] << " ";
        cout << "\n";
        cout << "Eigenvectors (rows) which are known in advance, Rij = \n";
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++)
            cout << Rt[i][j] << " ";
          cout << "\n";
        }
      }

      assert(Similar(evals, evals_known, n));
      for (int i = 0; i < n; i++)
        assert(Similar(R[i], Rt[i], n));

    } //for (int i_test = 0; i_test < n_tests_per_matrix; i++)

  } //for(int imat = 0; imat < n_matrices; imat++) {

  Dealloc2D(&_M, &M);
  Dealloc2D(&_R, &R);
  Dealloc2D(&_Rt, &Rt);
  Dealloc2D(&_D, &D);
  Dealloc2D(&_tmp, &tmp);
  delete [] evals;
  delete [] evals_known;
}


int main(int argc, char **argv) {
  int n_size = 2;
  int n_matr = 1;
  int n_tests = 1;
  Scalar erange = 1.0;

  if (argc <= 1) {
    cerr <<
      "Error: This program requires at least 1 arguments.\n"
      "Usage: n_size [n_matr n_tests erange]\n"
      "           (The 2nd and 3rd arguments are optional.)\n"
      "       n_size  = the size of the matrices\n"
      "       n_matr  = the number of randomly generated matrices to test\n"
      "       n_tests = the number of times the eigenvalues and eigenvectors\n"
      "                 are calculated for each matrix.  By default this is 1\n"
      "                 Increase it to at least 20 if you plan to use this\n"
      "                 program for benchmarking (speed testing), because the time\n"
      "                 needed for generating a random matrix is not negligible.\n"
      "       erange  = the range of eigenvalue magnitudes (log base 10), 1 by default\n"
         << endl;
    return 1;
  }

  n_size = std::stoi(argv[1]);
  if (argc >= 2)
    n_matr = std::stoi(argv[2]);
  if (argc >= 3)
    erange = std::stof(argv[3]);
  if (argc >= 4)
    n_tests = std::stoi(argv[4]);

  TestJacobi(n_size, n_matr, erange, n_tests);

  return EXIT_SUCCESS;
}
