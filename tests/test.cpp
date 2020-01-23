#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <random>
#include "matrix_alloc.hpp"
#include "jacobi.hpp"

using std::cout;
using std::cerr;
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
template<typename Scalar, typename Vector>
inline static bool SimilarVec(Vector a, Vector b, int n, Scalar eps=1.0e-06) {
  for (int i = 0; i < n; i++)
    if (not Similar(a[i], b[i], eps))
      return false;
  return true;
}

//Sort the rows of a matrix "evec" by the numbers contained in "eval"
//(This is a simple O(n^2) sorting method, but O(n^2) is a lower bound anyway.)
//This is the same as the Jacobi::SortRows(), but that function is private.
template<typename Scalar, typename Vector, typename Matrix>
void
SortRows(Vector eval, Matrix evec, int n, bool sort_absv=false)
{
  for (int i = 0; i < n; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
      if (sort_absv) { //sort by absolute value?
        if (std::abs(eval[j]) > std::abs(eval[i_max]))
          i_max = j;
      }
      else if (eval[j] > eval[i_max])
        i_max = j;
    }
    std::swap(eval[i], eval[i_max]); // sort "eval"
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]); // sort "evec"
  }
}



/// @brief
/// Generate a random orthogonal n x n matrix
template<typename Scalar, typename Matrix>
void GenRandOrth(Matrix R,
                 int n,
                 std::default_random_engine &generator)
{
  std::normal_distribution<double> gaussian_distribution(0,1);
  std::vector<Scalar> v(n);

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
        Scalar v_dot_Ri = 0.0;
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
    Scalar r_inv = 1.0 / std::sqrt(rsq);
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
                Scalar eval_magnitude_range=2.0, //<! range of eigevalues
                int n_tests_per_matrix=50) //<! repeat test for benchmarking?
{
  cout << endl << "-- Diagonalization test (real symmetric)  --" << endl;

  // Convert from base 10 to base e and (divide by 2)
  eval_magnitude_range *= std::log(10) * 0.5;

  Jacobi<double, double*, double**> eigen_calc(n);

  double *evals_known = new double[n];
  double *evals = new double[n];

  // construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> random_real01;

  double **M, **R, **Rt, **D, **tmp;
  double *_M, *_R, *_Rt, *_D, *_tmp;
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
      double rand_erange = (random_real01(generator)-0.5); // a number between [-0.5,0.5]
      rand_erange *= eval_magnitude_range; // scale this by eval_magnitude_range
      evals_known[i] = std::exp(rand_erange);
      // also consider both positive and negative eigenvalues:
      if (random_real01(generator) < 0.5)
        evals_known[i] = -evals_known[i];
      D[i][i] = evals_known[i];
    }

    // Now randomly generate the R and Rt matrices (ie. the eigenvectors):
    GenRandOrth<double, double**>(R, n, generator);
    // Rt is the transpose of R
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        Rt[i][j] = R[j][i];

    // Construct the test matrix, M, where M = Rt * D * R
    mmult(Rt, D, tmp, n);
    mmult(tmp, R, M, n);

    // Sort the matrix evals and eigenvector rows
    SortRows<double, double*, double**> (evals_known, R, n);

    if (n_matrices == 1) {
      cout << "Eigenvalues (after sorting):\n";
      for (int i = 0; i < n; i++)
        cout << evals_known[i] << " ";
      cout << "\n";
      cout << "Eigenvectors (rows) which are known in advance:\n";
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
          cout << R[i][j] << " ";
        cout << "\n";
      }
      cout << "  (The eigenvectors calculated by Jacobi::Diagonalize() should match these.)\n";
    }

    for (int i_test = 0; i_test < n_tests_per_matrix; i_test++) {

      // Now, calculate the eigenvalues and eigenvectors
      M[0][0]=0.65; //<-- CONTINUEHERE: for debugging only.  remove later!
      M[0][1]=0.53; //<-- CONTINUEHERE: for debugging only.  remove later!
      M[0][2]=0.11; //<-- CONTINUEHERE: for debugging only.  remove later!
      M[1][1]=0.5;  //<-- CONTINUEHERE: for debugging only.  remove later!
      M[1][2]=0.15; //<-- CONTINUEHERE: for debugging only.  remove later!
      M[2][2]=0.1;  //<-- CONTINUEHERE: for debugging only.  remove later!
      int n_sweeps = eigen_calc.Diagonalize(M, evals, Rt);

      if ((n_matrices == 1) && (i_test == 0)) {
        cout <<"Jacobi::Diagonalize() ran for "<<n_sweeps<<" iters (sweeps).\n";
        cout << "Eigenvalues calculated by Jacobi::Diagonalize()\n";
        for (int i = 0; i < n; i++)
          cout << evals[i] << " ";
        cout << "\n";
        cout << "Eigenvectors calculated by Jacobi::Diagonalize()\n";
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++)
            cout << Rt[i][j] << " ";
          cout << "\n";
        }
      }

      assert(SimilarVec(evals, evals_known, n, 1.0e-06));
      for (int i = 0; i < n; i++)
        assert(SimilarVec(R[i], Rt[i], n, 1.0e-06));

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
  double erange = 2.0;

  if (argc <= 1) {
    cerr <<
      "Error: This program requires at least 1 argument.\n"
      "Usage: n_size [n_matr n_tests erange]\n"
      "           (The remaining 3 arguments are optional.)\n"
      "       n_size  = the size of the matrices\n"
      "       n_matr  = the number of randomly generated matrices to test\n"
      "       n_tests = the number of times the eigenvalues and eigenvectors\n"
      "                 are calculated for each matrix.  By default this is 1\n"
      "                 (Increase this to at least 20 if you plan to use this\n"
      "                 program for benchmarking (speed testing), because the time\n"
      "                 needed for generating a random matrix is not negligible.)\n"
      "       erange  = the range of eigenvalue magnitudes (log base 10)\n"
      "#                A value of 2 (default) produces eigenvalues in the range from"
      "                 1.0 up to 10.0 (This is a 2-orders-of-magnitude difference.)"
      "                 A value of 1 produces eigenvalues from 1/sqrt(10) to sqrt(10)."
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
