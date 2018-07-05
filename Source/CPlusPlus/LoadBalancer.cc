#include "LoadBalancer.h"
#include "PMatrixMemoryPool.h"
#include "PSMatrix.h"
#include "Permutation.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "LoadBalancer_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
void LoadBalancer::PermuteMatrix(const Matrix_ps &mat_in, Matrix_ps &mat_out,
                                 const Permutation &permutation,
                                 PMatrixMemoryPool &memorypool) {
  PermuteMatrix_wrp(mat_in.ih_this, mat_out.ih_this, permutation.ih_this,
                    memorypool.ih_this);
}
void LoadBalancer::UndoPermuteMatrix(const Matrix_ps &mat_in,
                                     Matrix_ps &mat_out,
                                     const Permutation &permutation,
                                     PMatrixMemoryPool &memorypool) {
  UndoPermuteMatrix_wrp(mat_in.ih_this, mat_out.ih_this, permutation.ih_this,
                        memorypool.ih_this);
}
} // namespace NTPoly
