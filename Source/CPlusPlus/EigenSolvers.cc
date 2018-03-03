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
    DistributedSparseMatrix &eigenvectors,
    DistributedSparseMatrix &eigenvalues,
    const IterativeSolverParameters &solver_parameters) {
  EigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors), GetIH(eigenvalues),
                         GetIH(solver_parameters));
}
} // namespace NTPoly
