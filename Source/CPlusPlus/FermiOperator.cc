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
                                    const Matrix_ps &Overlap, int nel,
                                    Matrix_ps &Density, double inv_temp,
                                    double &energy_value_out,
                                    double &chemical_potential_out,
                                    const SolverParameters &solver_parameters) {
  ComputeDenseFOE_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
                      &inv_temp, &energy_value_out, &chemical_potential_out,
                      GetIH(solver_parameters));
}

} // namespace NTPoly
