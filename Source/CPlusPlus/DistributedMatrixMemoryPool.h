#ifndef DISTRIBUTEDMATRIXMEMORYPOOL_h
#define DISTRIBUTEDMATRIXMEMORYPOOL_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class DistributedSparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A memory pool datatype that can be reused for matrix matrix multiplication
//! this is to prevent excessive alloc/dealloc.
class DistributedMatrixMemoryPool {
public:
  // Standard constructor.
  DistributedMatrixMemoryPool();
  // Standard destructor.
  ~DistributedMatrixMemoryPool();

private:
  //! Pointer to underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  DistributedMatrixMemoryPool(const DistributedMatrixMemoryPool &);
  //! Assignment operator, locked.
  DistributedMatrixMemoryPool &operator=(const DistributedMatrixMemoryPool &);

private:
  friend class DistributedSparseMatrix;
  friend class LoadBalancer;
};
} // namespace NTPoly
#endif
