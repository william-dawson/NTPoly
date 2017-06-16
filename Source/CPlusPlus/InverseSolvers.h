#ifndef INVERSESOLVERS_h
#define INVERSESOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A class for inverting matrices.
class InverseSolvers : public SolverBase {
public:
  //! Compute the inverse of a matrix.
  //! An implementation of Hotelling's method.
  //!\param Overlap the matrix to invert.
  //!\param InverseMat = Overlap^-1.
  //!\param solver_parameters parameters for the solver
  static void Invert(const DistributedSparseMatrix &Overlap,
                     DistributedSparseMatrix &InverseMat,
                     const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
