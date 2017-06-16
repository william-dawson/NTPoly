#ifndef ROOTSOLVERS_h
#define ROOTSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing General Matrix Roots.
class RootSolvers : public SolverBase {
public:
  //! Compute a general matrix root.
  //!\param InputMat input matrix.
  //!\param OutputMat = InputMat^(1/root)
  //!\param root = which root
  //!\param solver_parameters parameters for the solver
  static void ComputeRoot(const DistributedSparseMatrix &InputMat,
                          DistributedSparseMatrix &OutputMat, int root,
                          const IterativeSolverParameters &solver_parameters);
  //! Compute the general matrix inverse root.
  //!\param InputMat input matrix.
  //!\param OutputMat = InputMat^(-1/root)
  //!\param root = which root
  //!\param solver_parameters parameters for the solver
  static void
  ComputeInverseRoot(const DistributedSparseMatrix &InputMat,
                     DistributedSparseMatrix &OutputMat, int root,
                     const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
