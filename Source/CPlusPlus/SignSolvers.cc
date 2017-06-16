#include "SignSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SignSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void SignSolvers::Compute(const DistributedSparseMatrix &mat1,
                          DistributedSparseMatrix &SignMat,
                          const IterativeSolverParameters &solver_parameters) {
  SignFunction_wrp(GetIH(mat1), GetIH(SignMat), GetIH(solver_parameters));
}
} // namespace NTPoly
