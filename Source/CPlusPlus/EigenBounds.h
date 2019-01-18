#ifndef EIGENBOUNDS_h
#define EIGENBOUNDS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing eigen bounds matrices.
class EigenBounds : public SolverBase {
public:
  //! Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  //! Uses Gershgorin's theorem.
  //!\param matrix the matrix to compute the min/max of.
  //!\param min_ger_eig a lower bound on the eigenspectrum.
  //!\param max_ger_eig an uppder bound on the eigenspectrum.
  static void GershgorinBounds(const Matrix_ps &matrix, double *min_ger_eig,
                               double *max_ger_eig);
  //! Compute a bounds on and maximum eigenvalue of a matrix.
  //! Uses The Power Method.
  //!\param matrix the matrix to compute the min/max of.
  //!\param max_power_eig an upper bound on the eigenspectrum.
  //!\param solver_parameters parameters for the solver
  static void PowerBounds(const Matrix_ps &matrix, double *max_power_eig,
                          const SolverParameters &solver_parameters);
  //! Compute interior eigenvalues of a matrix.
  //!\param matrix the matrix to compute the eigenvectors of.
  //!\param density The density matrix that splits the spectrum of this matrix.
  //!\param nel The number of electrons.
  //!\param nvals The number of values to compute. Negative if they should be
  //!below
  //! the gap, positive if above.
  //!\param vecs the output matrix of eigenvectors.
  //!\param solver_parameters parameters for the solver
  static void InteriorEigenvalues(const Matrix_ps &matrix,
                                  const Matrix_ps &density, int nel, int nvals,
                                  Matrix_ps &vecs,
                                  const SolverParameters &solver_parameters);
  //! Compute K largest eigenvalues with subspace iteration.
  //!\param matrix the matrix to compute the eigenvectors of.
  //!\param vecs the output matrix of eigenvectors.
  //!\param k the number of vectors to compute.
  //!\param solver_parameters parameters for the solver
  static void SubspaceIteration(const Matrix_ps &matrix, Matrix_ps &vecs, int k,
                                const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
