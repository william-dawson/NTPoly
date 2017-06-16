#ifndef LOADBALANCER_h
#define LOADBALANCER_h

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class DistributedMatrixMemoryPool;
class DistributedSparseMatrix;
class Permutation;

////////////////////////////////////////////////////////////////////////////////
//! A data structure for storing permutations.
class LoadBalancer {
public:
  //! Apply a permutation to a matrix.
  //\param[in] mat_in matrix to permute.
  //\param[out] mat_out permuted matrix.
  //\param[in] permutation to apply.
  //\param[inout] memorypool_in memory pool to use. Optional.
  static void PermuteMatrix(const DistributedSparseMatrix &mat_in,
                            DistributedSparseMatrix &mat_out,
                            const Permutation &permutation,
                            DistributedMatrixMemoryPool &memorypool);
  //! Undo a permutation applied to a matrix.
  //\param[in] mat_in matrix to undo permutation of.
  //\param[out] mat_out unpermuted matrix.
  //\param[in] permutation to remove.
  //\param[inout] memorypool_in memory pool to use. Optional.
  static void UndoPermuteMatrix(const DistributedSparseMatrix &mat_in,
                                DistributedSparseMatrix &mat_out,
                                const Permutation &permutation,
                                DistributedMatrixMemoryPool &memorypool);
};
} // namespace NTPoly
#endif
