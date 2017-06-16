#ifndef TRIGONOMETRY_h
#define TRIGONOMETRY_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class FixedSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Module For Computing Trigonometric functions of a Matrix.
class TrigonometrySolvers : public SolverBase {
public:
  //! Compute the sine of a matrix.
  //!\param InputMat matrix to compute the sine of.
  //!\param OutputMat = sin(InputMat)
  //!\param solver_parameters parameters for the solver
  static void Sine(const DistributedSparseMatrix &InputMat,
                   DistributedSparseMatrix &OutputMat,
                   const FixedSolverParameters &solver_parameters);
  //! Compute the cosine of a matrix.
  //!\param InputMat matrix to compute the cosine of.
  //!\param OutputMat = cos(InputMat)
  //!\param solver_parameters parameters for the solver
  static void Cosine(const DistributedSparseMatrix &InputMat,
                     DistributedSparseMatrix &OutputMat,
                     const FixedSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
