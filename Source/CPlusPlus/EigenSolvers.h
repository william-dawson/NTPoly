#ifndef EIGENSOLVERS_h
#define EIGENSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class DistributedSparseMatrix;
class FixedSolverParameters;
class IterativeSolverParameters;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing the eigen and singular value decomposition.
class EigenSolvers : public SolverBase {
public:
  //! Compute the eigenvalues and eigenvectors of a matrix
  static void ReferenceEigenDecomposition(
      const DistributedSparseMatrix &matrix,
      DistributedSparseMatrix &eigenvectors,
      DistributedSparseMatrix &eigenvalues,
      const FixedSolverParameters &solver_parameters);
  //! Compute the eigenvalues and eigenvectors of a matrix
  static void SplittingEigenDecomposition(
      const DistributedSparseMatrix &matrix,
      DistributedSparseMatrix &eigenvectors,
      DistributedSparseMatrix &eigenvalues, int num_values,
      const IterativeSolverParameters &solver_parameters);
  //! Compute the singularvalues and singularvectors of a matrix
  static void
  SingularValueDecompostion(const DistributedSparseMatrix &matrix,
                            DistributedSparseMatrix &leftvectors,
                            DistributedSparseMatrix &rightvectors,
                            DistributedSparseMatrix &singularvalues,
                            const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
