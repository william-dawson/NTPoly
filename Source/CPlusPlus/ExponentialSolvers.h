#ifndef EXPONENTIAL_h
#define EXPONENTIAL_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A Module For Computing General Matrix Exponentials.
class ExponentialSolvers : public SolverBase {
public:
  //! Compute the matrix exponential.
  //!\param InputMat matrix to compute the exponential of.
  //!\param OutputMat = exp(InputMat)
  //!\param solver_parameters parameters for the solver
  static void ComputeExponential(const Matrix_ps &InputMat,
                                 Matrix_ps &OutputMat,
                                 const SolverParameters &solver_parameters);
  //! Compute the matrix exponential using the Pade method.
  //!\param InputMat matrix to compute the exponential of.
  //!\param OutputMat = exp(InputMat)
  //!\param solver_parameters parameters for the solver
  static void ComputeExponentialPade(const Matrix_ps &InputMat,
                                     Matrix_ps &OutputMat,
                                     const SolverParameters &solver_parameters);
  //! Compute the matrix logarithm.
  //!\param InputMat matrix to compute the exponential of.
  //!\param OutputMat = log(InputMat)
  //!\param solver_parameters parameters for the solver
  static void ComputeLogarithm(const Matrix_ps &InputMat, Matrix_ps &OutputMat,
                               const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
