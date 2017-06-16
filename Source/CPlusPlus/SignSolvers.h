#ifndef SIGNSOLVERS_h
#define SIGNSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A Class For Computing The Matrix Sign Function.
class SignSolvers : public SolverBase {
public:
  //! Compute the matrix sign function.
  //!\param Mat1 input matrix.
  //!\param SignMat = Sign(Mat1)
  //!\param solver_parameters parameters for the solver
  static void Compute(const DistributedSparseMatrix &Mat1,
                      DistributedSparseMatrix &SignMat,
                      const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
