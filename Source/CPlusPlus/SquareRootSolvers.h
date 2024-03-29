#ifndef SQUAREROOT_h
#define SQUAREROOT_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing The Square Root of a Matrix.
class SquareRootSolvers : public SolverBase {
public:
  //! Compute the square root of a matrix.
  //!\param InputMat matrix to compute the inversesquareroot of.
  //!\param OutputMat = InputMat^1/2.
  //!\param solver_parameters parameters for the solver
  static void SquareRoot(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                         const SolverParameters &solver_parameters);
  //! Compute the square root of a matrix (dense version).
  //!\param InputMat matrix to compute the inversesquareroot of.
  //!\param OutputMat = InputMat^1/2.
  //!\param solver_parameters parameters for the solver
  static void DenseSquareRoot(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                              const SolverParameters &solver_parameters);
  //! Compute the inverse square root of a matrix.
  //!\param InputMat matrix to compute the inversesquareroot of.
  //!\param OutputMat = InputMat^-1/2.
  //!\param solver_parameters parameters for the solver
  static void InverseSquareRoot(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                                const SolverParameters &solver_parameters);
  //! Compute the inverse square root of a matrix (dense version).
  //!\param InputMat matrix to compute the inversesquareroot of.
  //!\param OutputMat = InputMat^-1/2.
  //!\param solver_parameters parameters for the solver
  static void DenseInverseSquareRoot(const Matrix_ps &InputMat,
                                     Matrix_ps &OutputMat,
                                     const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
