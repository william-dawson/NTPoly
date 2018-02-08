#include "SignSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SignSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void SignSolvers::ComputeSign(
    const DistributedSparseMatrix &mat1, DistributedSparseMatrix &SignMat,
    const IterativeSolverParameters &solver_parameters) {
  SignFunction_wrp(GetIH(mat1), GetIH(SignMat), GetIH(solver_parameters));
}
void SignSolvers::ComputePolarDecomposition(
    const DistributedSparseMatrix &mat1, DistributedSparseMatrix &Umat,
    DistributedSparseMatrix &Hmat,
    const IterativeSolverParameters &solver_parameters) {
  PolarDecomposition_wrp(GetIH(mat1), GetIH(Umat), GetIH(Hmat),
                         GetIH(solver_parameters));
}
} // namespace NTPoly
