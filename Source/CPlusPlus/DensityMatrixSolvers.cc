#include "DensityMatrixSolvers.h"
#include "PSMatrix.h"
#include "SolverParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "DensityMatrixSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::PM(const Matrix_ps &Hamiltonian,
                              const Matrix_ps &Overlap, int nel,
                              Matrix_ps &Density, double &energy_value_out,
                              double &chemical_potential_out,
                              const SolverParameters &solver_parameters) {
  PM_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
         &energy_value_out, &chemical_potential_out, GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::TRS2(const Matrix_ps &Hamiltonian,
                                const Matrix_ps &Overlap, int nel,
                                Matrix_ps &Density, double &energy_value_out,
                                double &chemical_potential_out,
                                const SolverParameters &solver_parameters) {
  TRS2_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::TRS4(const Matrix_ps &Hamiltonian,
                                const Matrix_ps &Overlap, int nel,
                                Matrix_ps &Density, double &energy_value_out,
                                double &chemical_potential_out,
                                const SolverParameters &solver_parameters) {
  TRS4_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::HPCP(const Matrix_ps &Hamiltonian,
                                const Matrix_ps &Overlap, int nel,
                                Matrix_ps &Density, double &energy_value_out,
                                double &chemical_potential_out,
                                const SolverParameters &solver_parameters) {
  HPCP_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
// void DensityMatrixSolvers::HPCPPlus(
//     const Matrix_ps &Hamiltonian,
//     const Matrix_ps &Overlap, int nel,
//     Matrix_ps &Density, double &energy_value_out, double
//     &chemical_potential_out, const SolverParameters &solver_parameters) {
//   HPCPPlus_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
//                &energy_value_out, &energy_value_out, &chemical_potential_out,
//                GetIH(solver_parameters));
// }

} // namespace NTPoly
