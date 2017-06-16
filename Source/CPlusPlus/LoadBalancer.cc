#include "DistributedBlockedSparseMatrix.h"
#include "DistributedMatrixMemoryPool.h"
#include "LoadBalancer.h"
#include "Permutation.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "LoadBalancer_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
void LoadBalancer::PermuteMatrix(const DistributedSparseMatrix &mat_in,
                                 DistributedSparseMatrix &mat_out,
                                 const Permutation &permutation,
                                 DistributedMatrixMemoryPool &memorypool) {
  PermuteMatrix_wrp(mat_in.ih_this, mat_out.ih_this, permutation.ih_this,
                    memorypool.ih_this);
}
void LoadBalancer::UndoPermuteMatrix(const DistributedSparseMatrix &mat_in,
                                     DistributedSparseMatrix &mat_out,
                                     const Permutation &permutation,
                                     DistributedMatrixMemoryPool &memorypool) {
  UndoPermuteMatrix_wrp(mat_in.ih_this, mat_out.ih_this, permutation.ih_this,
                        memorypool.ih_this);
}
} // namespace NTPoly
