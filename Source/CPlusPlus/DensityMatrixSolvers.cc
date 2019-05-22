#include "DensityMatrixSolvers.h"
#include "PSMatrix.h"
#include "SolverParameters.h"
using std::vector;
#include <iostream>

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
void DensityMatrixSolvers::ScaleAndFold(
    const Matrix_ps &Hamiltonian, const Matrix_ps &Overlap, int nel,
    Matrix_ps &Density, double homo, double lumo, double &energy_value_out,
    const SolverParameters &solver_parameters) {
  ScaleAndFold_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
                   &homo, &lumo, &energy_value_out, GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::DACPurify(
    const Matrix_ps &Hamiltonian, const Matrix_ps &Overlap,
    vector<int> gap_list, Matrix_ps &Density, double &energy_value_out,
    double &chemical_potential_out, const SolverParameters &solver_parameters) {
  int len;
  len = (int)gap_list.size();
  DACPurify_wrp(GetIH(Hamiltonian), GetIH(Overlap), &gap_list[0], &len,
                GetIH(Density), &energy_value_out, &chemical_potential_out,
                GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::EnergyDensityMatrix(const Matrix_ps &Hamiltonian,
                                               const Matrix_ps &Density,
                                               Matrix_ps &EnergyDensity,
                                               double threshold) {
  EnergyDensityMatrix_wrp(GetIH(Hamiltonian), GetIH(Density),
                          GetIH(EnergyDensity), &threshold);
}

} // namespace NTPoly
