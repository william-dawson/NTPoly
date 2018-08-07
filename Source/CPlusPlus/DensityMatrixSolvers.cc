#include "DensityMatrixSolvers.h"
#include "DistributedSparseMatrix.h"
#include "IterativeSolversParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "DensityMatrixSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::PM(
    const DistributedSparseMatrix &Hamiltonian,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &Density, double &energy_value_out,
    double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  PM_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
         &energy_value_out, &chemical_potential_out, GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::TRS2(
    const DistributedSparseMatrix &Hamiltonian,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &Density, double &energy_value_out,
    double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  TRS2_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::TRS4(
    const DistributedSparseMatrix &Hamiltonian,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &Density, double &energy_value_out,
    double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  TRS4_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void DensityMatrixSolvers::HPCP(
    const DistributedSparseMatrix &Hamiltonian,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &Density, double &energy_value_out,
    double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  HPCP_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
           &energy_value_out, &chemical_potential_out,
           GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
// void DensityMatrixSolvers::HPCPPlus(
//     const DistributedSparseMatrix &Hamiltonian,
//     const DistributedSparseMatrix &Overlap, int nel,
//     DistributedSparseMatrix &Density, double &energy_value_out,double
//     &chemical_potential_out, const IterativeSolverParameters
//     &solver_parameters) {
//   HPCPPlus_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
//                &energy_value_out, &chemical_potential_out,
//                GetIH(solver_parameters));
// }

} // namespace NTPoly
