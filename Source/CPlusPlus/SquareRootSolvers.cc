#include "SquareRootSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SquareRootSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void SquareRootSolvers::SquareRoot(const Matrix_ps &Input, Matrix_ps &Output,
                                   const SolverParameters &solver_parameters) {
  SquareRoot_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void SquareRootSolvers::InverseSquareRoot(
    const Matrix_ps &Input, Matrix_ps &Output,
    const SolverParameters &solver_parameters) {
  InverseSquareRoot_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
} // namespace NTPoly
