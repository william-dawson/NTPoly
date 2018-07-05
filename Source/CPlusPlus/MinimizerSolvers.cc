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
    Matrix_ps &Density, double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  ConjugateGradient_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel,
                        GetIH(Density), &chemical_potential_out,
                        GetIH(solver_parameters));
}
} // namespace NTPoly
