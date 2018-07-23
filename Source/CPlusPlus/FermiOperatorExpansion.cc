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
    Matrix_ps &Density, int degree, double chemical_potential,
    const SolverParameters &solver_parameters) {
  ComputeFOEFixed_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
                      &degree, &chemical_potential, GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void FermiOperatorExpansion::Compute(
    const Matrix_ps &Hamiltonian, const Matrix_ps &Overlap, int nel,
    Matrix_ps &Density, int degree, const SolverParameters &solver_parameters) {
  ComputeFOESearch_wrp(GetIH(Hamiltonian), GetIH(Overlap), &nel, GetIH(Density),
                       &degree, GetIH(solver_parameters));
}

////////////////////////////////////////////////////////////////////////////////
void FermiOperatorExpansion::ComputeEigenvalues(
    const Matrix_ps &InputMat, const Matrix_ps &Overlap,
    const Matrix_ps &Eigenvalues, int degree, int nvals,
    const SolverParameters &solver_parameters) {
  FOEEigenvalues_wrp(GetIH(InputMat), GetIH(Overlap), GetIH(Eigenvalues),
                     &degree, &nvals, GetIH(solver_parameters));
}
} // namespace NTPoly
