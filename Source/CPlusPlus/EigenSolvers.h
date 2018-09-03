#ifndef EIGENSOLVERS_h
#define EIGENSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
class SolverParameters;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing the eigen and singular value decomposition.
class EigenSolvers : public SolverBase {
public:
  //! Compute the eigenvalues and eigenvectors of a matrix
  static void
  ReferenceEigenDecomposition(const Matrix_ps &matrix, Matrix_ps &eigenvectors,
                              Matrix_ps &eigenvalues,
                              const SolverParameters &solver_parameters);
  //! Compute the eigenvalues and eigenvectors of a matrix
  static void
  SplittingEigenDecomposition(const Matrix_ps &matrix, Matrix_ps &eigenvectors,
                              Matrix_ps &eigenvalues, int num_values,
                              const SolverParameters &solver_parameters);
  //! Compute the singularvalues and singularvectors of a matrix
  static void
  SingularValueDecompostion(const Matrix_ps &matrix, Matrix_ps &leftvectors,
                            Matrix_ps &rightvectors, Matrix_ps &singularvalues,
                            const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
