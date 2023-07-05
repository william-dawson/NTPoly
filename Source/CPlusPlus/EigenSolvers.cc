#include "EigenSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "EigenSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void EigenSolvers::EigenDecomposition(
    const Matrix_ps &matrix, Matrix_ps &eigenvalues, int nvals,
    Matrix_ps &eigenvectors, const SolverParameters &solver_parameters) {
  EigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvalues), &nvals,
                         GetIH(eigenvectors), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void EigenSolvers::EigenValues(const Matrix_ps &matrix, Matrix_ps &eigenvalues,
                               int nvals,
                               const SolverParameters &solver_parameters) {
  EigenDecomposition_novec_wrp(GetIH(matrix), GetIH(eigenvalues), &nvals,
                               GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void EigenSolvers::SingularValueDecomposition(
    const Matrix_ps &matrix, Matrix_ps &leftvectors, Matrix_ps &rightvectors,
    Matrix_ps &singularvalues, const SolverParameters &solver_parameters) {
  SingularValueDecompostion_wrp(GetIH(matrix), GetIH(leftvectors),
                                GetIH(rightvectors), GetIH(singularvalues),
                                GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void EigenSolvers::EstimateGap(const Matrix_ps &H, const Matrix_ps &K,
                               double chemical_potential, double *gap,
                               const SolverParameters &solver_parameters) {
  EstimateGap_wrp(GetIH(H), GetIH(K), &chemical_potential, gap,
                  GetIH(solver_parameters));
}
} // namespace NTPoly
