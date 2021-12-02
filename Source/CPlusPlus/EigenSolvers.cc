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
    const Matrix_ps &matrix, Matrix_ps &eigenvectors, Matrix_ps &eigenvalues,
    const SolverParameters &solver_parameters) {
  EigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors), GetIH(eigenvalues),
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
} // namespace NTPoly
