#include "EigenSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "EigenSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void EigenSolvers::ReferenceEigenDecomposition(
    const Matrix_ps &matrix, Matrix_ps &eigenvectors, Matrix_ps &eigenvalues,
    const FixedSolverParameters &solver_parameters) {
  ReferenceEigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors),
                                  GetIH(eigenvalues), GetIH(solver_parameters));
}
void EigenSolvers::SplittingEigenDecomposition(
    const Matrix_ps &matrix, Matrix_ps &eigenvectors, Matrix_ps &eigenvalues,
    int num_values, const IterativeSolverParameters &solver_parameters) {
  SplittingEigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors),
                                  GetIH(eigenvalues), &num_values,
                                  GetIH(solver_parameters));
}
void EigenSolvers::SingularValueDecompostion(
    const Matrix_ps &matrix, Matrix_ps &leftvectors, Matrix_ps &rightvectors,
    Matrix_ps &singularvalues,
    const IterativeSolverParameters &solver_parameters) {
  SingularValueDecompostion_wrp(GetIH(matrix), GetIH(leftvectors),
                                GetIH(rightvectors), GetIH(singularvalues),
                                GetIH(solver_parameters));
}
} // namespace NTPoly
