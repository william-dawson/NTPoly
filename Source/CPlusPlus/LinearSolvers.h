#ifndef LINEARSOLVERS_h
#define LINEARSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class FixedSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A class for solving matrix equations.
class LinearSolvers : public SolverBase {
public:
  //! Solve the matrix equation AX = B using conjugate gradient method.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param XMat the solved for matrix X.
  //!\param BMat the right hand side.
  //!\param solver_parameters parameters for the solver
  static void CGSolver(const DistributedSparseMatrix &AMat,
                       DistributedSparseMatrix &XMat,
                       const DistributedSparseMatrix &BMat,
                       const IterativeSolverParameters &solver_parameters);
  //! Compute The Cholesky Decomposition of a Symmetric Positive Definite
  //! matrix.
  //!\param AMat the matrix A, must be symmetric, positive definite.
  //!\param LMat the matrix computed.
  //!\param solver_parameters parameters for the solver
  static void
  CholeskyDecomposition(const DistributedSparseMatrix &AMat,
                        DistributedSparseMatrix &LMat,
                        const FixedSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
