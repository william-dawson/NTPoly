#include "DistributedSparseMatrix.h"
#include "GeometryOptimization.h"
#include "IterativeSolversParameters.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "GeometryOptimization_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void GeometryOptimization::PurificationExtrapolate(
    const DistributedSparseMatrix &PreviousDensity,
    const DistributedSparseMatrix &Overlap, int nel,
    DistributedSparseMatrix &NewDensity,
    const IterativeSolverParameters &solver_parameters) {
  PurificationExtrapolate_wrp(GetIH(PreviousDensity), GetIH(Overlap), &nel,
                              GetIH(NewDensity), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void GeometryOptimization::LowdinExtrapolate(
    const DistributedSparseMatrix &PreviousDensity,
    const DistributedSparseMatrix &OldOverlap,
    const DistributedSparseMatrix &NewOverlap,
    DistributedSparseMatrix &NewDensity,
    const IterativeSolverParameters &solver_parameters) {
  LowdinExtrapolate_wrp(GetIH(PreviousDensity), GetIH(OldOverlap),
                        GetIH(NewOverlap), GetIH(NewDensity),
                        GetIH(solver_parameters));
}
} // namespace NTPoly
