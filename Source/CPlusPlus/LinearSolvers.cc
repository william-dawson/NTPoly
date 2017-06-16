#include "LinearSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "LinearSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void LinearSolvers::CGSolver(
    const DistributedSparseMatrix &AMat, DistributedSparseMatrix &XMat,
    const DistributedSparseMatrix &BMat,
    const IterativeSolverParameters &solver_parameters) {
  CGSolver_wrp(GetIH(AMat), GetIH(XMat), GetIH(BMat), GetIH(solver_parameters));
}
} // namespace NTPoly
