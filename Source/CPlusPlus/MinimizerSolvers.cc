#include "MinimizerSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "MinimizerSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void MinimizerSolvers::ConjugateGradient(
    const Matrix_ps &Hamiltonian, const Matrix_ps &Overlap, int nel,
    Matrix_ps &Density, double &energy_value_out,
    double &chemical_potential_out, const SolverParameters &solver_parameters) {
  ConjugateGradient_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel,
                        GetIH(Density), &energy_value_out,
                        &chemical_potential_out, GetIH(solver_parameters));
}
} // namespace NTPoly
