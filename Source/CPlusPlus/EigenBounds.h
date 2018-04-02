#ifndef EIGENBOUNDS_h
#define EIGENBOUNDS_h

#include "SolverBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing eigen bounds matrices.
class EigenBounds : public SolverBase {
public:
  //! Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  //! Uses Gershgorin's theorem.
  //!\param matrix the matrix to compute the min/max of.
  //!\param min_ger_eig a lower bound on the eigenspectrum.
  //!\param max_ger_eig an uppder bound on the eigenspectrum.
  static void GershgorinBounds(const DistributedSparseMatrix &matrix,
                               double *min_ger_eig, double *max_ger_eig);
  //! Compute a bounds on and maximum eigenvalue of a matrix.
  //! Uses The Power Method.
  //!\param matrix the matrix to compute the min/max of.
  //!\param max_power_eig an upper bound on the eigenspectrum.
  //!\param solver_parameters parameters for the solver
  static void PowerBounds(const DistributedSparseMatrix &matrix,
                          double *max_power_eig,
                          const IterativeSolverParameters &solver_parameters);

  static void
  EigenDecomposition(const DistributedSparseMatrix &matrix,
                     DistributedSparseMatrix &eigenvectors,
                     const FixedSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
