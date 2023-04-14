#include "FermiOperator.h"
#include "PSMatrix.h"
#include "SolverParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "FermiOperator_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void FermiOperator::ComputeDenseFOE(const Matrix_ps &Hamiltonian,
                                    const Matrix_ps &Overlap, double trace,
                                    Matrix_ps &Density, double inv_temp,
                                    double &energy_value_out,
                                    double &chemical_potential_out,
                                    const SolverParameters &solver_parameters) {
  ComputeDenseFOE_wrp(GetIH(Hamiltonian), GetIH(Overlap), &trace,
                      GetIH(Density), &inv_temp, &energy_value_out,
                      &chemical_potential_out, GetIH(solver_parameters));
}
void FermiOperator::WOM_GC(const Matrix_ps &Hamiltonian,
                           const Matrix_ps &Overlap,
                           Matrix_ps &Density, double chemical_potential, 
                           double inv_temp, double &energy_value_out,
                           const SolverParameters &solver_parameters) {
  WOM_GC_wrp(GetIH(Hamiltonian), GetIH(Overlap),
             GetIH(Density), &chemical_potential, &inv_temp, 
             &energy_value_out, GetIH(solver_parameters));
}

} // namespace NTPoly
