#ifndef GEOMETRYOPTIMIZATIONSOLVERS_ch
#define GEOMETRYOPTIMIZATIONSOLVERS_ch

void PurificationExtrapolate_wrp(const int *ih_PreviousDensity,
                                 const int *Overlap, const int *nel,
                                 int *ih_NewDensity,
                                 const int *ih_solver_parameters);
#endif
