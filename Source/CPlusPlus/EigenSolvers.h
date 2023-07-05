#ifndef EIGENSOLVERS_h
#define EIGENSOLVERS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing eigen decompositions matrices.
class EigenSolvers : public SolverBase {
public:
  //! Compute the eigendecomposition of a matrix.
  //! Uses a dense routine.
  //!\param matrix the matrix to decompose.
  //!\param eigenvalues the eigenvalues of a matrix.
  //!\param eigenvectors the eigenvectors of a matrix.
  //!\param nvals the number of values to compute.
  //!\param solver_parameters parameters for computing.
  static void EigenDecomposition(const Matrix_ps &matrix,
                                 Matrix_ps &eigenvalues, int nvals,
                                 Matrix_ps &eigenvectors,
                                 const SolverParameters &solver_parameters);
  //! Compute the eigenvalues of a matrix.
  //! Uses a dense routine.
  //!\param matrix the matrix to decompose.
  //!\param eigenvalues the eigenvalues of a matrix.
  //!\param nvals the number of values to compute.
  //!\param solver_parameters parameters for computing.
  static void EigenValues(const Matrix_ps &matrix, Matrix_ps &eigenvalues,
                          int nvals, const SolverParameters &solver_parameters);
  //! Compute the singular value decomposition of a matrix.
  //! Uses a dense routine.
  //!\param matrix the matrix to decompose.
  //!\param leftvectors the left singular vectors.
  //!\param rightvectors the right singular vectors
  //!\param singularvalues a diagonal matrix containing the singular values.
  //!\param solver_parameters parameters for computing.
  static void
  SingularValueDecomposition(const Matrix_ps &matrix, Matrix_ps &leftvectors,
                             Matrix_ps &rightvectors, Matrix_ps &singularvalues,
                             const SolverParameters &solver_parameters);
  //! Estimate the HOMO-LUMO gap of a matrix.
  //!\param H matrix to compute the gap of.
  //!\param ISQ inverse square root of the overlap matrix.
  //!\param K density matrix.
  //!\param chemical_potential (estimated from purification).
  //!\param gap computed gap value.
  //!\param solver_parameters parameters for computing.
  static void EstimateGap(const Matrix_ps &H, const Matrix_ps &K,
                          double chemical_potential, double* gap,
                          const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
