#ifndef EXPONENTIAL_h
#define EXPONENTIAL_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class FixedSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Module For Computing General Matrix Exponentials.
class ExponentialSolvers : public SolverBase {
public:
  //! Compute the matrix exponential.
  //!\param InputMat matrix to compute the exponential of.
  //!\param OutputMat = exp(InputMat)
  //!\param solver_parameters parameters for the solver
  static void
  ComputeExponential(const DistributedSparseMatrix &InputMat,
                     DistributedSparseMatrix &OutputMat,
                     const FixedSolverParameters &solver_parameters);
  //! Compute the matrix logarithm.
  //!\param InputMat matrix to compute the exponential of.
  //!\param OutputMat = log(InputMat)
  //!\param solver_parameters parameters for the solver
  static void ComputeLogarithm(const DistributedSparseMatrix &InputMat,
                               DistributedSparseMatrix &OutputMat,
                               const FixedSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
