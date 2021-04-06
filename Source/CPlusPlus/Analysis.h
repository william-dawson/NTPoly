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
};
} // namespace NTPoly
#endif
