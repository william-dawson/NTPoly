#include "SquareRootSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SquareRootSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void SquareRootSolvers::SquareRoot(
    const DistributedSparseMatrix &Input, DistributedSparseMatrix &Output,
    const IterativeSolverParameters &solver_parameters) {
  SquareRoot_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void SquareRootSolvers::InverseSquareRoot(
    const DistributedSparseMatrix &Input, DistributedSparseMatrix &Output,
    const IterativeSolverParameters &solver_parameters) {
  InverseSquareRoot_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
} // namespace NTPoly
