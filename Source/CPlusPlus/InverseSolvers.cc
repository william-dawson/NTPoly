#include "InverseSolvers.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "InverseSolvers_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
void InverseSolvers::CholeskyInvert(const Matrix_ps &Matrix,
                                    Matrix_ps &InverseMat,
                                    const SolverParameters &solver_parameters) {
  CholeskyInvert_wrp(GetIH(Matrix), GetIH(InverseMat),
                     GetIH(solver_parameters));
}
void InverseSolvers::Invert(const Matrix_ps &Matrix, Matrix_ps &InverseMat,
                            const SolverParameters &solver_parameters) {
  Invert_wrp(GetIH(Matrix), GetIH(InverseMat), GetIH(solver_parameters));
}
////////////////////////////////////////////////////////////////////////////////
void InverseSolvers::PseudoInverse(const Matrix_ps &Matrix,
                                   Matrix_ps &InverseMat,
                                   const SolverParameters &solver_parameters) {
  PseudoInverse_wrp(GetIH(Matrix), GetIH(InverseMat), GetIH(solver_parameters));
}
} // namespace NTPoly
