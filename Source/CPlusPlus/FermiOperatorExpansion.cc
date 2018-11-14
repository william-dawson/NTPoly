#include "FermiOperatorExpansion.h"
#include "PSMatrix.h"
#include "SolverParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "FermiOperatorExpansion_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void FermiOperatorExpansion::Compute(
    const Matrix_ps &Hamiltonian, const Matrix_ps &Overlap, int nel,
    Matrix_ps &Density, int degree, double &energy_value_out,
    double &chemical_potential_out, const SolverParameters &solver_parameters) {
  ComputeFOE_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
                 &degree, &energy_value_out, &chemical_potential_out,
                 GetIH(solver_parameters));
}

} // namespace NTPoly
