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
    const DistributedSparseMatrix &Hamiltonian,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &Density, double &chemical_potential_out,
    const IterativeSolverParameters &solver_parameters) {
  ConjugateGradient_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel,
                        GetIH(Density), &chemical_potential_out,
                        GetIH(solver_parameters));
}
} // namespace NTPoly
