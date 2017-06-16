#include "TrigonometrySolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "TrigonometrySolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void TrigonometrySolvers::Sine(const DistributedSparseMatrix &Input,
                               DistributedSparseMatrix &Output,
                               const FixedSolverParameters &solver_parameters) {
  Sine_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void TrigonometrySolvers::Cosine(
    const DistributedSparseMatrix &Input, DistributedSparseMatrix &Output,
    const FixedSolverParameters &solver_parameters) {
  Cosine_wrp(GetIH(Input), GetIH(Output), GetIH(solver_parameters));
}
} // namespace NTPoly
