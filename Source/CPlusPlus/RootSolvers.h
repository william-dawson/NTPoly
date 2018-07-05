#ifndef ROOTSOLVERS_h
#define ROOTSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing General Matrix Roots.
class RootSolvers : public SolverBase {
public:
  //! Compute a general matrix root.
  //!\param InputMat input matrix.
  //!\param OutputMat = InputMat^(1/root)
  //!\param root = which root
  //!\param solver_parameters parameters for the solver
  static void ComputeRoot(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                          int root, const SolverParameters &solver_parameters);
  //! Compute the general matrix inverse root.
  //!\param InputMat input matrix.
  //!\param OutputMat = InputMat^(-1/root)
  //!\param root = which root
  //!\param solver_parameters parameters for the solver
  static void ComputeInverseRoot(const Matrix_ps &InputMat,
                                 Matrix_ps &OutputMat, int root,
                                 const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
