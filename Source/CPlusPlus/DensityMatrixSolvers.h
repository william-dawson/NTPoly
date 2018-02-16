#ifndef DENSITYMATRIXSOLVERS_h
#define DENSITYMATRIXSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class IterativeSolverParameters;
class DistributedSparseMatrix;
//! A Class For Solving Chemistry Systems Based On Sparse Matrices.
class DensityMatrixSolvers : public SolverBase {
public:
  //! Compute the density matrix from a Hamiltonian using the TRS2 method.
  //! Based on the TRS2 algorithm presented in:
  //! Niklasson, Anders MN. "Expansion algorithm for the density matrix."
  //! Physical Review B 66, no. 15 (2002): 155115.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void TRS2(const DistributedSparseMatrix &Hamiltonian,
                   const DistributedSparseMatrix &InverseSquareRoot, int nel,
                   DistributedSparseMatrix &Density,
                   double &chemical_potential_out,
                   const IterativeSolverParameters &solver_parameters);
  //! Compute the density matrix from a Hamiltonian using the TRS4 method.
  //! Based on the TRS4 algorithm presented in:
  //! Niklasson, Anders MN. "Expansion algorithm for the density matrix."
  //! Physical Review B 66, no. 15 (2002): 155115.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void TRS4(const DistributedSparseMatrix &Hamiltonian,
                   const DistributedSparseMatrix &InverseSquareRoot, int nel,
                   DistributedSparseMatrix &Density,
                   double &chemical_potential_out,
                   const IterativeSolverParameters &solver_parameters);
  //! Compute the density matrix from a Hamiltonian using the TDB method.
  //! Based on the algorithm presented in:
  //! Truflandier, Lionel A., Rivo M. Dianzinga, and David R. Bowler.
  //! "Communication: Generalized canonical purification for density matrix
  //! minimization." The Journal of chemical physics 144, no. 9 (2016): 091102
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void HPCP(const DistributedSparseMatrix &Hamiltonian,
                   const DistributedSparseMatrix &InverseSquareRoot, int nel,
                   DistributedSparseMatrix &Density,
                   double &chemical_potential_out,
                   const IterativeSolverParameters &solver_parameters);

  //! Create a new guess at the Density Matrix after updating the geometry.
  //! Based on the purification algorithm in \cite niklasson2010trace .
  //!\param PreviousDensity to extrapolate from.
  //!\param Overlap the overlap matrix of the new geometry.
  //!\param nel the number of electrons.
  //!\param NewDensity the extrapolated density.
  //!\param solver_parameters parameters for the solver
  static void
  ExtrapolateGeometry(const DistributedSparseMatrix &PreviousDensity,
                      const DistributedSparseMatrix &Overlap, int nel,
                      DistributedSparseMatrix &NewDensity,
                      const IterativeSolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
