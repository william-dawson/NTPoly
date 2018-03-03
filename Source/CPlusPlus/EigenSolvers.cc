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
    const IterativeSolverParameters &solver_parameters) {
  EigenDecomposition_wrp(GetIH(matrix), GetIH(eigenvectors),
                         GetIH(solver_parameters));
}
} // namespace NTPoly
