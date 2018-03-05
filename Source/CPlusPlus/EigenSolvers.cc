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
    const DistributedSparseMatrix &matrix,
    DistributedSparseMatrix &eigenvectors, DistributedSparseMatrix &eigenvalues,
    int num_values, const IterativeSolverParameters &solver_parameters) {
  EigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors), GetIH(eigenvalues),
                         &num_values, GetIH(solver_parameters));
}
void EigenSolvers::SingularValueDecompostion(
    const DistributedSparseMatrix &matrix, DistributedSparseMatrix &leftvectors,
    DistributedSparseMatrix &rightvectors,
    DistributedSparseMatrix &singularvalues,
    const IterativeSolverParameters &solver_parameters) {
  SingularValueDecompostion_wrp(GetIH(matrix), GetIH(leftvectors),
                                GetIH(rightvectors), GetIH(singularvalues),
                                GetIH(solver_parameters));
}
} // namespace NTPoly
