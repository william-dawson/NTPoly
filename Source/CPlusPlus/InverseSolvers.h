#ifndef INVERSESOLVERS_h
#define INVERSESOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for inverting matrices.
class InverseSolvers : public SolverBase {
public:
  //! Compute the inverse of a matrix using the Cholesky decomposition.
  //!\param Matrix the matrix to invert.
  //!\param InverseMat = Matrix^-1.
  //!\param solver_parameters parameters for the solver
  static void CholeskyInvert(const Matrix_ps &Matrix, Matrix_ps &InverseMat,
                     const SolverParameters &solver_parameters);
  //! Compute the inverse of a matrix.
  //! An implementation of Hotelling's method.
  //!\param Matrix the matrix to invert.
  //!\param InverseMat = Matrix^-1.
  //!\param solver_parameters parameters for the solver
  static void Invert(const Matrix_ps &Matrix, Matrix_ps &InverseMat,
                     const SolverParameters &solver_parameters);
  //! Compute the pseudoinverse of a matrix.
  //! An implementation of Hotelling's method, with a different convergence
  //! criteria.
  //!\param Matrix the matrix to invert.
  //!\param InverseMat = Matrix^-1.
  //!\param solver_parameters parameters for the solver
  static void PseudoInverse(const Matrix_ps &Matrix, Matrix_ps &InverseMat,
                            const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
