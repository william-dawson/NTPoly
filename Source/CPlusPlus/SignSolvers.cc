#include "SignSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SignSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void SignSolvers::ComputeSign(const Matrix_ps &mat1, Matrix_ps &SignMat,
                              const SolverParameters &solver_parameters) {
  SignFunction_wrp(GetIH(mat1), GetIH(SignMat), GetIH(solver_parameters));
}
void SignSolvers::ComputePolarDecomposition(
    const Matrix_ps &mat1, Matrix_ps &Umat, Matrix_ps &Hmat,
    const SolverParameters &solver_parameters) {
  PolarDecomposition_wrp(GetIH(mat1), GetIH(Umat), GetIH(Hmat),
                         GetIH(solver_parameters));
}
} // namespace NTPoly
