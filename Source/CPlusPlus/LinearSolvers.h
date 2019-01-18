#ifndef LINEARSOLVERS_h
#define LINEARSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for solving matrix equations.
class LinearSolvers : public SolverBase {
public:
  //! Solve the matrix equation AX = B using conjugate gradient method.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param XMat the solved for matrix X.
  //!\param BMat the right hand side.
  //!\param solver_parameters parameters for the solver
  static void CGSolver(const Matrix_ps &AMat, Matrix_ps &XMat,
                       const Matrix_ps &BMat,
                       const SolverParameters &solver_parameters);
  //! Compute The Cholesky Decomposition of a Symmetric Positive Definite
  //! matrix.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param LMat the matrix computed.
  //!\param solver_parameters parameters for the solver
  static void CholeskyDecomposition(const Matrix_ps &AMat, Matrix_ps &LMat,
                                    const SolverParameters &solver_parameters);
  //! Compute The Cholesky Decomposition of a Symmetric Positive Semi-Definite
  //! matrix.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param LMat the matrix computed.
  //!\param rank the target rank
  //!\param solver_parameters parameters for the solver
  static void
  PivotedCholeskyDecomposition(const Matrix_ps &AMat, Matrix_ps &LMat, int rank,
                               const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
