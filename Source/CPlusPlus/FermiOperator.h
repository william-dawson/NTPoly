#ifndef FERMIOPERATOREXPANSION_h
#define FERMIOPERATOREXPANSION_h

#include "SolverBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SolverParameters;
class Matrix_ps;
//! A Class For Solving Chemistry Systems Using the Fermi Operator Expansion.
class FermiOperator : public SolverBase {
public:
  //! Compute the density matrix using a dense solver.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param trace of the density matrix (usually the number of electrons).
  //!\param Density the density matrix computed by this routine.
  //!\param inv_temp the inverse temperature.
  //!\param energy_value_out the energy of the system.
  //!\param chemical_potential_out the chemical potential calculated.
  //!\param solver_parameters parameters for the solver
  static void ComputeDenseFOE(const Matrix_ps &Hamiltonian,
                              const Matrix_ps &InverseSquareRoot, double trace,
                              Matrix_ps &Density, double inv_temp,
                              double &energy_value_out,
                              double &chemical_potential_out,
                              const SolverParameters &solver_parameters);
  //! Compute the density matrix using wave operator minimization.
  //! Fixed chemical potential version.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param Density the density matrix computed by this routine.
  //!\param chemical_potential the chemical potential calculated.
  //!\param inv_temp the inverse temperature.
  //!\param energy_value_out the energy of the system.
  //!\param solver_parameters parameters for the solver
  static void WOM_GC(const Matrix_ps &Hamiltonian,
                     const Matrix_ps &InverseSquareRoot, Matrix_ps &Density,
                     double chemical_potential, double inv_temp,
                     double &energy_value_out,
                     const SolverParameters &solver_parameters);
  //! Compute the density matrix using wave operator minimization.
  //! Fixed number of electrons version.
  //!\param Hamiltonian the matrix to compute the corresponding density from.
  //!\param InverseSquareRoot of the overlap matrix.
  //!\param Density the density matrix computed by this routine.
  //!\param trace the target trace of the matrix.
  //!\param inv_temp the inverse temperature.
  //!\param energy_value_out the energy of the system.
  //!\param solver_parameters parameters for the solver
  static void WOM_C(const Matrix_ps &Hamiltonian,
                    const Matrix_ps &InverseSquareRoot, Matrix_ps &Density,
                    double trace, double inv_temp, double &energy_value_out,
                    const SolverParameters &solver_parameters);
};
} // namespace NTPoly
#endif
