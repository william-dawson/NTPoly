#include "EigenBounds.h"
using namespace NTPoly;
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "EigenBounds_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void EigenBounds::GershgorinBounds(
    const DistributedSparseMatrix& matrix, double * min_value,
    double * max_value) {
  GershgorinBounds_wrp(GetIH(matrix), min_value, max_value);
}
void EigenBounds::PowerBounds(
    const DistributedSparseMatrix &matrix, double * max_value,
    const IterativeSolverParameters &solver_parameters) {
  PowerBounds_wrp(GetIH(matrix), max_value, GetIH(solver_parameters));
}
} // namespace NTPoly
