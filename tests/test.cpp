#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <random>
#include "matrix_alloc.hpp"

using std::cout;
using std::endl;
using std::setprecision;
using namespace matrix_alloc;

template<typename T>
using vector = std::vector<T>;

template<typename T>
using complex = std::complex<T>;

/// @brief
/// Generate a random orthogonal n x n matrix
template<typename Scalar, typename Matrix>
void GenRandOrth(int n,
	         Matrix R)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0,1); //(Gaussian distribution)

  for (int i = 0; i < n; i++) {
    // Generate a vector in a random direction
    // (This works because we are using a normal (a.k.a. Gaussian) distribution)
    for (int j = 0; j < n; j++)
    v[j] = distribution(generator);

    //Now subtract from v, the projection of v onto the first i-1 rows of R.
    //This will produce a vector which is orthogonal to these i-1 row-vectors.
    for (int k = 0; k < i; k++) {
      Scalar v_dot_Ri = 0.0
      for (int j = 0; j < n; j++)
        v_dot_Ri += v[j] * R[k][j]; // = <v , R[i]>
      for (int j = 0; j < n; j++)
        v[j] -= v_dot_Ri * R[i][j] // = v - <V,R[i]> R[i]
    }
    // Now normalize what remains
    rsq = 0.0;
    for (int j = 0; j < n; j++)
      rsq += v[j]*v[j];
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




/// @brief  Convert a 3x3 rotation matrix (M) to a quaternion (q)
/// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

template<typename Scalar>
static inline void Matrix2Quaternion(const Scalar M[3][3], Scalar q[4])
{
  Scalar S;
  Scalar qw, qx, qy, qz;
  Scalar tr = Trace3(M);  // = M[0][0] + M[1][1] + M[2][2];
  if (tr > 0) {
    S = std::sqrt(tr+1.0) * 2;                        // S=4*qw 
    qw = 0.25 * S;
    qx = (M[2][1] - M[1][2]) / S;
    qy = (M[0][2] - M[2][0]) / S;
    qz = (M[1][0] - M[0][1]) / S;
  }
  else if ((M[0][0] > M[1][1]) and (M[0][0] > M[2][2])) {
    S = std::sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2;   // S=4*qx 
    qw = (M[2][1] - M[1][2]) / S;
    qx = 0.25 * S;
    qy = (M[0][1] + M[1][0]) / S;
    qz = (M[0][2] + M[2][0]) / S;
  }
  else if (M[1][1] > M[2][2]) {
    S = std::sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2;   // S=4*qy
    qw = (M[0][2] - M[2][0]) / S;
    qx = (M[0][1] + M[1][0]) / S;
    qy = 0.25 * S;
    qz = (M[1][2] + M[2][1]) / S;
  }
  else {
    S = std::sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2;   // S=4*qz
    qw = (M[1][0] - M[0][1]) / S;
    qx = (M[0][2] + M[2][0]) / S;
    qy = (M[1][2] + M[2][1]) / S;
    qz = 0.25 * S;
  }
  q[0] = qw;
  q[1] = qx;
  q[2] = qy;
  q[3] = qz;

} //Matrix2Quaternion(const Scalar M[3][3], Scalar q[4])






/// @brief  Convert a quaternion (q) to a 3x3 rotation matrix (M)
/// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm

template<typename Scalar>
static inline void Quaternion2Matrix(Scalar const *q, Scalar **M)
{
  // Alternate convention for quaternion->rotation conversion
  // This could be the convention where q[0] is real instead of q[3] (?I forgot)
  //M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  //M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  //M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3]);
  //M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  //M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  //M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  //M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  //M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
  //M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);

  // Quaternion->rotation conversion used here:
  // This could be the convention where q[3] is real instead of q[0] (?I forgot)
  M[0][0] =  1.0 - 2*(q[2]*q[2]) - 2*(q[3]*q[3]);
  M[1][1] =  1.0 - 2*(q[1]*q[1]) - 2*(q[3]*q[3]);
  M[2][2] =  1.0 - 2*(q[1]*q[1]) - 2*(q[2]*q[2]);
  M[0][1] = 2*(q[1]*q[2] - q[3]*q[0]);
  M[1][0] = 2*(q[1]*q[2] + q[3]*q[0]);
  M[1][2] = 2*(q[2]*q[3] - q[1]*q[0]);
  M[2][1] = 2*(q[2]*q[3] + q[1]*q[0]);
  M[0][2] = 2*(q[1]*q[3] + q[2]*q[0]);
  M[2][0] = 2*(q[1]*q[3] - q[2]*q[0]);
}


/// @brief  Convert a 3-component Shoemake coordinate list to a quaternion (q)
///         Shoemake, Graphics Gems III (1992) pp. 124-132
template <typename Scalar>
static inline void Shoemake2Quaternion(const Scalar sm[3], Scalar q[4])
{
  const Scalar M_2PI = 6.283185307179586;
  Scalar X0 = sm[0];
  Scalar X1 = sm[1];
  Scalar X2 = sm[2];
  Scalar theta1 = M_2PI * X1;
  Scalar theta2 = M_2PI * X2;
  Scalar r1 = std::sqrt(1.0-X0);
  Scalar r2 = std::sqrt(X0);
  Scalar s1 = std::sin(theta1);
  Scalar c1 = std::cos(theta1);
  Scalar s2 = std::sin(theta2);
  Scalar c2 = std::cos(theta2);

  // Alternative quaternion convention, where q[3] is real instead of q[0]
  q[0] = s1*r1;
  q[1] = c1*r1;
  q[2] = s2*r2;
  q[3] = c2*r2;

  // Alternative quaternion convention, where q[0] is real instead of q[3]
  //q[3] = s1*r1;
  //q[0] = c1*r1;
  //q[1] = s2*r2;
  //q[2] = c2*r2;
}



/// @brief  Convert a quaternion (q) to a 3-component Shoemake coordinate list
///         Shoemake, Graphics Gems III (1992) pp. 124-132
template <typename Scalar>
static inline void Quaternion2Shoemake(const Scalar q[4], Scalar sm[3])
{
  const Scalar M_2PI = 6.283185307179586;

  // Alternative quaternion convention, where q[3] is real instead of q[0]
  Scalar r1 = std::sqrt(q[0]*q[0] + q[1]*q[1]);
  Scalar r2 = std::sqrt(q[2]*q[2] + q[3]*q[3]);
  Scalar X0 = r2*r2;
  Scalar theta1 = 0.0;
  if (r1 > 0)
    theta1 = std::atan2(q[0], q[1]);
  Scalar theta2 = 0.0;
  if (r2 > 0)
    theta2 = std::atan2(q[2], q[3]);

  // Alternative quaternion convention, where q[0] is real instead of q[3]
  //Scalar r1 = std::sqrt(q[3]*q[3] + q[0]*q[0]);
  //Scalar r2 = std::sqrt(q[1]*q[1] + q[2]*q[2]);
  //Scalar X0 = r2;
  //Scalar theta1 = 0.0;
  //if (r1 > 0)
  //  theta1 = std::atan2(q[3], q[0]);
  //Scalar theta2 = 0.0;
  //if (r2 > 0)
  //  theta2 = std::atan2(q[1], q[2]);

  Scalar X1 = theta1 / M_2PI;
  Scalar X2 = theta2 / M_2PI;
  sm[0] = X0;
  sm[1] = X1;
  sm[2] = X2;
}



template<typename Scalar>
void Quaternion2MatrixALT(Scalar const *q, Scalar **M) {
  // convert a quaternion (q) to a 3x3 rotation matrix (M)
  M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3]);
  M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3]);
  M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
  M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);
}

template<typename Scalar>
void RotMatAXYZ(Scalar **dest,
                Scalar angle,
                Scalar axis_x,
                Scalar axis_y,
                Scalar axis_z)
{

  Scalar r = std::sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);

  Scalar X = 1.0;
  Scalar Y = 0.0;
  Scalar Z = 0.0;
  if (r > 0.0) {
    //check for non-sensical input
    X = axis_x / r;
    Y = axis_y / r;
    Z = axis_z / r;
  }
  else
    angle = 0.0;

  // angle *= math.pi/180.0 # "angle" is assumed to be in degrees
  //    on second thought, let the caller worry about angle units.
  Scalar c = std::cos(angle);
  Scalar s = std::sin(angle);

  dest[0][0] = X * X * (1 - c) + c;
  dest[1][1] = Y * Y * (1 - c) + c;
  dest[2][2] = Z * Z * (1 - c) + c;

  dest[0][1] = X * Y * (1 - c) - Z * s;
  dest[0][2] = X * Z * (1 - c) + Y * s;

  dest[1][0] = Y * X * (1 - c) + Z * s;
  dest[2][0] = Z * X * (1 - c) - Y * s;

  dest[1][2] = Y * Z * (1 - c) - X * s;
  dest[2][1] = Z * Y * (1 - c) + X * s;

  //   formula from these sources:
  // http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
  // also check
  // http://www.manpagez.com/man/3/glRotate/
  //   some pdb test commands:
  // from lttree_matrixstack import *
  // r = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  // RotMatAXYZ(r, 90.0, 0.0, 0.0, 1.0)
}



/// @brief  Convert 3 Shoemake coordinates to a 3x3 rotation matrix (M)
template <typename Scalar>
static inline void Shoemake2Matrix(const Scalar sm[3], Scalar M[3][3]) {
  Scalar q[4];
  Shoemake2Quaternion(sm, q);
  Quaternion2Matrix(q, M);
}


/// @brief  Convert 3 Shoemake coordinates to a 3x3 rotation matrix (M)
template <typename Scalar>
static inline void Matrix2Shoemake(const Scalar M[3][3], Scalar sm[3]) {
  Scalar q[4];
  Matrix2Quaternion(M, q);
  Quaternion2Shoemake(q, sm);
}



void test1() {
  cout << endl << "-- Diagonalization test (real symmetric)  --" << endl;
  
  const int n = 3;
  //double matrix[n][n] = { {2.0, 1.0, 1.0},
  //                        {1.0, 2.0, 1.0},
  //                        {1.0, 1.0, 2.0} };
  ///* Its eigenvalues are {4, 1, 1} */

  double **M, **R, **Rt, **D, **tmp;
  double *_M, *_R, *_Rt, *_D, *_tmp;
  Alloc2D(3, 3, &_M, &M);
  Alloc2D(3, 3, &_R, &R);
  Alloc2D(3, 3, &_Rt, &Rt);
  Alloc2D(3, 3, &_D, &D);
  Alloc2D(3, 3, &_tmp, &tmp);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      D[i][j] = 0.0;
  D[0][0] = 100.0;
  D[1][1] = 101.00001;
  D[2][2] = 101.0;
  RotMatAXYZ(R, -2.7, 0.1, 0.9, -0.7);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Rt[i][j] = R[j][i];
  
  MatProduct3(D, Rt, tmp);
  MatProduct3(R, tmp, M);

  cout << "Eigenvectors (columns) which are known in advance, Rij = \n";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << R[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "  (The eigenvector should match one of these columns.)\n";
  cout << "\n";

  PEigenCalculator<double> pe(3, true);
  double evec[3];
  double eval = pe.PrincipalEigen(M, evec);
  cout << " --- Results using PEigenCalculator ---\n";
  cout << "Eigen value: " << setprecision(16) << eval << endl;
  cout << "Eigen vector:";
  for(int i = 0;i < n;i++) {
    cout << evec[i] << " ";
  }
  cout << endl;
  
  return EXIT_SUCCESS;
}
