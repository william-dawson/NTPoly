#ifndef ANALYSIS_h
#define ANALYSIS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for solving matrix equations.
class Analysis : public SolverBase {
public:
  //! Compute The Cholesky Decomposition of a Symmetric Positive Semi-Definite
  //! matrix.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param LMat the matrix computed.
  //!\param rank the target rank
  //!\param solver_parameters parameters for the solver
  static void
  PivotedCholeskyDecomposition(const Matrix_ps &AMat, Matrix_ps &LMat, int rank,
                               const SolverParameters &solver_parameters);
  //!  When we want to only compute the first n eigenvalues of a matrix, this
  //!  routine will project out the higher eigenvalues.
  //!\param AMat The starting matrix.
  //!\param dim  The number of eigenvalues ot keep.
  //!\param RMat a dimxdim matrix with the same first n eigenvalues as A.
  //!\param solver_parameters parameters for the solver
  static void
  ReduceDimension(const Matrix_ps &AMat, int dim, Matrix_ps &RMat,
                  const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
