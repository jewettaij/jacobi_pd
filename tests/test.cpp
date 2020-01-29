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

// @brief  Are two numbers "similar"?
template<typename T>
inline static bool Similar(T a, T b, T eps=1.0e-06) {
  return std::abs(a - b) <= std::abs(eps)*std::sqrt(std::abs(a)*std::abs(b));
}

/// @brief  Are two vectors (containing n numbers) similar?
template<typename Scalar, typename Vector>
inline static bool SimilarVec(Vector a, Vector b, int n, Scalar eps=1.0e-06) {
  for (int i = 0; i < n; i++)
    if (not Similar(a[i], b[i], eps))
      return false;
  return true;
}

/// @brief  Are two vectors (or their reflections) similar?
template<typename Scalar, typename Vector>
inline static bool SimilarVecUnsigned(Vector a, Vector b, int n, Scalar eps=1.0e-06) {
  if (SimilarVec(a, b, n, eps))
    return true;
  else {
    for (int i = 0; i < n; i++)
      if (not Similar(a[i], -b[i], eps))
        return false;
    return true;
  }
}


/// @brief  Multiply two matrices A and B, store the result in C. (C = AB).

template<typename Scalar, typename Matrix, typename ConstMatrix>
void mmult(ConstMatrix A, //<! input array
           ConstMatrix B, //<! input array
           Matrix C,      //<! store result here
           int m,      //<! number of rows of A
           int n=0,    //<! optional: number of columns of B (=m by default)
           int K=0     //<! optional: number of columns of A = num rows of B (=m by default)
           )
{
  assert((C != A) && (C != B));
  if (n == 0) n = m; // if not specified, then assume the matrices are square
  if (K == 0) K = m; // if not specified, then assume the matrices are square

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      C[i][j] = 0.0;

  // perform matrix multiplication
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < K; k++)
        C[i][j] += A[i][k] * B[k][j];
}



/// @brief
///Sort the rows of a matrix "evec" by the numbers contained in "eval"
///(This is a simple O(n^2) sorting method, but O(n^2) is a lower bound anyway.)
///This is the same as the Jacobi::SortRows(), but that function is private.
template<typename Scalar, typename Vector, typename Matrix>
void
SortRows(Vector eval,
         Matrix evec,
         int n,
         bool sort_decreasing=true,
         bool sort_abs=false)
{
  for (int i = 0; i < n-1; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
      if (sort_decreasing) {
        if (sort_abs) { //sort by absolute value?
          if (std::abs(eval[j]) > std::abs(eval[i_max]))
            i_max = j;
        }
        else if (eval[j] > eval[i_max])
          i_max = j;
      }
      else {
        if (sort_abs) { //sort by absolute value?
          if (std::abs(eval[j]) < std::abs(eval[i_max]))
            i_max = j;
        }
        else if (eval[j] < eval[i_max])
          i_max = j;
      }
    }
    std::swap(eval[i], eval[i_max]); // sort "eval"
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]); // sort "evec"
  }
}



/// @brief  Generate a random orthonormal n x n matrix

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
      //(They are already normalized and orthogonal to each other.)
      for (int k = 0; k < i; k++) {
        Scalar v_dot_Rk = 0.0;
          for (int j = 0; j < n; j++)
            v_dot_Rk += v[j] * R[k][j];
        for (int j = 0; j < n; j++)
          v[j] -= v_dot_Rk * R[k][j];
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



/// @brief  Generate a random symmetric n x n matrix, M.
/// This function generates random numbers for the eigenvalues ("evals_known")
/// as well as the eigenvectors ("evects_known"), and uses them to generate M.
/// The "eval_magnitude_range" argument specifies the the base-10 logarithm
/// of the range of eigenvalues desired.  The "n_degeneracy" argument specifies
/// the number of repeated eigenvalues desired (if any).
/// @returns  This function does not return a value.  However after it is
///           invoked, the M matrix will be filled with random numbers.
///           Additionally, the "evals" and "evects" arguments will contain
///           the eigenvalues and eigenvectors (one eigenvector per row)
///           of the matrix.  Later, they can be compared with the eigenvalues
///           and eigenvectors calculated by Jacobi::Diagonalize()

template <typename Scalar, typename Vector, typename Matrix>
void GenRandSymm(Matrix M, //<! store the matrix here
                 int n, //<! matrix size
                 Scalar eval_magnitude_range=2.0,//<! range of wanted eigevalues
                 int n_degeneracy=1,//<!number of repeated eigevalues(1disables)
                 std::default_random_engine &generator,//<! makes random numbers
                 Vector evals,  //<! store the eigenvalues of here
                 Matrix evects  //<! store the eigenvectors here
                 )
{
  Scalar  **D, **tmp;
  Alloc2D(n, n, &D);
  Alloc2D(n, n, &tmp);

  // D is a diagonal matrix whose diagonal elements are the eigenvalues
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      D[i][j] = 0.0;

  // Randomly generate the eigenvalues
  for (int i = 0; i < n; i++) {
    // generate some random eigenvalues
    // Use a "log-uniform distribution" (a.k.a. "reciprocal distribution")
    // (This is a way to specify numbers with a precise range of magnitudes.)
    Scalar rand_erange = (random_real01(generator)-0.5); // a number between [-0.5,0.5]
    rand_erange *= eval_magnitude_range; // scale this by eval_magnitude_range
    evals[i] = std::exp(rand_erange*std::log(10.0));
    // also consider both positive and negative eigenvalues:
    if (random_real01(generator) < 0.5)
      evals[i] = -evals[i];
    D[i][i] = evals[i];
  }

  // Now randomly generate the R and Rt matrices (ie. the eigenvectors):
  GenRandOrth<Scalar, Matrix>(evects, n, generator);

  // Construct the test matrix, M, where M = Rt * D * R
  mmult(evects, D, tmp, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evects[i][j] = evects[j][i]; //transpose "evects"
  mmult(tmp, evects, M, n);  //at this point M = Rt * D * R (where "R"="evects")
  Dealloc2D(&D);
  Dealloc2D(&tmp);
} // GenRandSymm()



template <typename Scalar>
void TestJacobi(int n, //<! matrix size
                int n_matrices=100, //<! number of matrices to test
                Scalar eval_magnitude_range=2.0, //<! range of eigevalues
                int n_tests_per_matrix=1, //<! repeat test for benchmarking?
                unsigned seed=0 //<! random seed (if 0 then use the clock)
                )
{
  cout << endl << "-- Diagonalization test (real symmetric)  --" << endl;

  // Convert from base 10 to base e and (divide by 2)
  eval_magnitude_range *= std::log(10) * 0.5;

  Jacobi<Scalar, Scalar*, Scalar**, Scalar const*const*> eigen_calc(n);

  // construct a random generator engine using a time-based seed:

  if (seed == 0) // if the caller did not specify a seed, use the system clock
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<Scalar> random_real01;

  Scalar **M, **evects, **evects_known;
  Alloc2D(n, n, &M);
  Alloc2D(n, n, &evects);
  Alloc2D(n, n, &evects_known);
  Scalar *evals = new Scalar[n];
  Scalar *evals_known = new Scalar[n];
  Scalar *test_evec = new Scalar[n];

  for(int imat = 0; imat < n_matrices; imat++) {

    // Create a randomly generated symmetric matrix.
    //This function generates random numbers for the eigenvalues ("evals_known")
    //as well as the eigenvectors ("evects_known"), and uses them to generate M.

    GenRandSymm<Scalar, Scalar*, Scalar**>(M,
                                           n,
                                           eval_magnitude_range,
                                           n_degeneracy,
                                           generator,
                                           evals_known,
                                           evects_known);

    // Sort the matrix evals and eigenvector rows
    SortRows<Scalar, Scalar*, Scalar**> (evals_known, evects_known, n);

    if (n_matrices == 1) {
      cout << "Eigenvalues (after sorting):\n";
      for (int i = 0; i < n; i++)
        cout << evals_known[i] << " ";
      cout << "\n";
      cout << "Eigenvectors (rows) which are known in advance:\n";
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
          cout << evects_known[i][j] << " ";
        cout << "\n";
      }
      cout << "  (The eigenvectors calculated by Jacobi::Diagonalize() should match these.)\n";
    }

    for (int i_test = 0; i_test < n_tests_per_matrix; i_test++) {

      // Now, calculate the eigenvalues and eigenvectors
      int n_sweeps = eigen_calc.Diagonalize(M, evals, evects);

      if ((n_matrices == 1) && (i_test == 0)) {
        cout <<"Jacobi::Diagonalize() ran for "<<n_sweeps<<" iters (sweeps).\n";
        cout << "Eigenvalues calculated by Jacobi::Diagonalize()\n";
        for (int i = 0; i < n; i++)
          cout << evals[i] << " ";
        cout << "\n";
        cout << "Eigenvectors (rows) calculated by Jacobi::Diagonalize()\n";
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++)
            cout << evects[i][j] << " ";
          cout << "\n";
        }
      }

      assert(SimilarVec(evals, evals_known, n, 1.0e-06));
      //Check each eigenvector
      for (int i = 0; i < n; i++) {
        for (int a = 0; a < n; a++) {
          test_evec[i] = 0.0;
          for (int b = 0; b < n; b++)
            test_evec[a] += M[a][b] * evec[i][b];
          assert(Similar(test_evec[a], eval[i] * evec[i][a], 1.0e-06));
        }
      }

    } //for (int i_test = 0; i_test < n_tests_per_matrix; i++)

  } //for(int imat = 0; imat < n_matrices; imat++) {

  Dealloc2D(&M);
  Dealloc2D(&evects_known);
  Dealloc2D(&evects);
  delete [] evals_known;
  delete [] evals;
  delete [] test_evec;

}


int main(int argc, char **argv) {
  int n_size = 2;
  int n_matr = 1;
  double erange = 2.0;
  int n_tests = 1;
  int n_degeneracy = 1;
  unsigned seed = 0;

  if (argc <= 1) {
    cerr <<
      "Error: This program requires at least 1 argument.\n"
      "Usage: n_size [n_matr n_tests erange n_degeneracy seed]\n"
      "       n_size  = the size of the matrices\n"
      "           (The remaining arguments are optional.)\n"
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
      "  n_degeneracy = the number of repeated eigenvalues (disabled by default)\n"
      "          seed = the seed used by the random number generator.\n"
      "                 (By default, the system clock is used.)\n"
         << endl;
    return 1;
  }

  n_size = std::stoi(argv[1]);
  if (argc > 2)
    n_matr = std::stoi(argv[2]);
  if (argc > 3)
    erange = std::stof(argv[3]);
  if (argc > 4)
    n_tests = std::stoi(argv[4]);
  if (argc > 6)
    n_degeneracy = std::stoi(argv[6]);
  if (argc > 7)
    seed = std::stoi(argv[7]);

  TestJacobi(n_size, n_matr, erange, n_tests, n_degeneracy, seed);

  cout << "test passed\n" << endl;
  return EXIT_SUCCESS;
}
