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
void LinearSolvers::CholeskyDecomposition(
    const DistributedSparseMatrix &AMat, DistributedSparseMatrix &LMat,
    const FixedSolverParameters &solver_parameters) {
  CholeskyDecomposition_wrp(GetIH(AMat), GetIH(LMat), GetIH(solver_parameters));
}
void LinearSolvers::PivotedCholeskyDecomposition(
    const DistributedSparseMatrix &AMat, DistributedSparseMatrix &LMat,
    int rank, const FixedSolverParameters &solver_parameters) {
  PivotedCholeskyDecomposition_wrp(GetIH(AMat), GetIH(LMat), &rank,
                                   GetIH(solver_parameters));
}
} // namespace NTPoly
