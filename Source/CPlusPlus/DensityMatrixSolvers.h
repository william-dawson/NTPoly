#ifndef DENSITYMATRIXSOLVERS_h
#define DENSITYMATRIXSOLVERS_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
//! A Class For Solving Chemistry Systems Based On Sparse Matrices.
class DensityMatrixSolvers : public SolverBase {
public:
  //! Compute the density matrix from a Hamiltonian using the PM method.
  //! Based on the PM algorithm presented in \cite palser1998canonical
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param energy_value_out the energy of the system (optional).
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void PM(const Matrix_ps &Hamiltonian,
                 const Matrix_ps &InverseSquareRoot, int nel,
                 Matrix_ps &Density, double &energy_value_out,
                 double &chemical_potential_out,
                 const SolverParameters &solver_parameters);
  //! Compute the density matrix from a Hamiltonian using the TRS2 method.
  //! Based on the TRS2 algorithm presented in: \cite niklasson2002.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param energy_value_out the energy of the system (optional).
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void TRS2(const Matrix_ps &Hamiltonian,
                   const Matrix_ps &InverseSquareRoot, int nel,
                   Matrix_ps &Density, double &energy_value_out,
                   double &chemical_potential_out,
                   const SolverParameters &solver_parameters);
  //! Compute the density matrix from a Hamiltonian using the TRS4 method.
  //! Based on the TRS4 algorithm presented in: \cite niklasson2002 .
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param energy_value_out the energy of the system (optional).
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void TRS4(const Matrix_ps &Hamiltonian,
                   const Matrix_ps &InverseSquareRoot, int nel,
                   Matrix_ps &Density, double &energy_value_out,
                   double &chemical_potential_out,
                   const SolverParameters &solver_parameters);
  //! Compute the density matrix from a Hamiltonian using the HPCP method.
  //! Based on the algorithm presented in: \cite truflandier2016communication
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param energy_value_out the energy of the system (optional).
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void HPCP(const Matrix_ps &Hamiltonian,
                   const Matrix_ps &InverseSquareRoot, int nel,
                   Matrix_ps &Density, double &energy_value_out,
                   double &chemical_potential_out,
                   const SolverParameters &solver_parameters);
  //! Compute the density matrix of a system using the eigendecomposition.
  //! Based on the algorithm presented in: \cite truflandier2016communication
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param nel the number of electrons.
  //!\param Density the density matrix computed by this routine.
  //!\param energy_value_out the energy of the system (optional).
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void DenseSolver(const Matrix_ps &Hamiltonian,
                          const Matrix_ps &InverseSquareRoot, int nel,
                          Matrix_ps &Density, double &energy_value_out,
                          double &chemical_potential_out,
                          const SolverParameters &solver_parameters);
  //! Compute the energy-weighted density matrix.
  //!\param Hamiltonian the matrix to compute from.
  //!\param Density the density matrix.
  //!\param EnergyDensity the energy-weighted density matrix to compute.
  //!\param threshold for flushing small values to zero.
  static void EnergyDensityMatrix(const Matrix_ps &Hamiltonian,
                                  const Matrix_ps &Density,
                                  Matrix_ps &EnergyDensity,
                                  double threshold = 0.0);
};
} // namespace NTPoly
#endif
