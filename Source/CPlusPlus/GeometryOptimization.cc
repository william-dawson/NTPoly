#include "GeometryOptimization.h"
#include "IterativeSolversParameters.h"
#include "PSMatrix.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "GeometryOptimization_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void GeometryOptimization::PurificationExtrapolate(
    const Matrix_ps &PreviousDensity, const Matrix_ps &Overlap, int nel,
    Matrix_ps &NewDensity, const IterativeSolverParameters &solver_parameters) {
  PurificationExtrapolate_wrp(GetIH(PreviousDensity), GetIH(Overlap), &nel,
                              GetIH(NewDensity), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void GeometryOptimization::LowdinExtrapolate(
    const Matrix_ps &PreviousDensity, const Matrix_ps &OldOverlap,
    const Matrix_ps &NewOverlap, Matrix_ps &NewDensity,
    const IterativeSolverParameters &solver_parameters) {
  LowdinExtrapolate_wrp(GetIH(PreviousDensity), GetIH(OldOverlap),
                        GetIH(NewOverlap), GetIH(NewDensity),
                        GetIH(solver_parameters));
}
} // namespace NTPoly
