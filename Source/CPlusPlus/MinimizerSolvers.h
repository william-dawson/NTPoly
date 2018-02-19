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
  //! Based on two papers. Based on two papers. The first by Scuseria
  //! \cite millam1997linear and then Challacombe
  //! \cite challacombe1999simplified
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
