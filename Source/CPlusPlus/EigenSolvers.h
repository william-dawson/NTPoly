#ifndef EIGENSOLVERS_h
#define EIGENSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing the eigen and singular value decomposition.
class EigenSolvers : public SolverBase {
public:
  //! Compute the eigenvalues and eigenvectors of a matrix
  static void
  EigenDecomposition(const DistributedSparseMatrix &matrix,
                     DistributedSparseMatrix &eigenvectors,
                     DistributedSparseMatrix &eigenvalues,
                     const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
