#ifndef MINIMIZERSOLVERS_h
#define MINIMIZERSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A class for computing the density matrix based on minimization methods.
class MinimizerSolvers : public SolverBase {
public:
  //! Compute the density matrix from a Hamiltonian using the CG method.
  //! Based on two papers. The first by Scuseria:
  //! Millam, John M., and Gustavo E. Scuseria. "Linear scaling conjugate
  //! gradient density matrix search as an alternative to diagonalization for
  //! first principles electronic structure calculations." The Journal of
  //! chemical physics 106, no. 13 (1997): 5569-5577. The second one by
  //! Chalacombe: Challacombe, Matt. "A simplified density matrix minimization
  //! for linear scaling self-consistent field theory." The Journal of chemical
  //! physics 110, no. 5 (1999): 2332-2342.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void
  ConjugateGradient(const DistributedSparseMatrix &Hamiltonian,
                    const DistributedSparseMatrix &InverseSquareRoot, int nel,
                    DistributedSparseMatrix &Density,
                    double &chemical_potential_out,
                    const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
