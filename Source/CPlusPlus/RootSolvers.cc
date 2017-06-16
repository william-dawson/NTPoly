#include "RootSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "RootSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void RootSolvers::ComputeRoot(
    const DistributedSparseMatrix &InputMat, DistributedSparseMatrix &OutputMat,
    int root, const IterativeSolverParameters &solver_parameters) {
  ComputeRoot_wrp(GetIH(InputMat), GetIH(OutputMat), &root,
                  GetIH(solver_parameters));
}
//////////////////////////////////////////////////////////////////////////////
void RootSolvers::ComputeInverseRoot(
    const DistributedSparseMatrix &InputMat, DistributedSparseMatrix &OutputMat,
    int root, const IterativeSolverParameters &solver_parameters) {
  ComputeInverseRoot_wrp(GetIH(InputMat), GetIH(OutputMat), &root,
                         GetIH(solver_parameters));
}
} // namespace NTPoly
