#ifndef LOADBALANCER_h
#define LOADBALANCER_h

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class PMatrixMemoryPool;
class Matrix_ps;
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
  static void PermuteMatrix(const Matrix_ps &mat_in, Matrix_ps &mat_out,
                            const Permutation &permutation,
                            PMatrixMemoryPool &memorypool);
  //! Undo a permutation applied to a matrix.
  //\param[in] mat_in matrix to undo permutation of.
  //\param[out] mat_out unpermuted matrix.
  //\param[in] permutation to remove.
  //\param[inout] memorypool_in memory pool to use. Optional.
  static void UndoPermuteMatrix(const Matrix_ps &mat_in, Matrix_ps &mat_out,
                                const Permutation &permutation,
                                PMatrixMemoryPool &memorypool);
};
} // namespace NTPoly
#endif
