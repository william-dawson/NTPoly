#ifndef DISTRIBUTEDMATRIXMEMORYPOOL_h
#define DISTRIBUTEDMATRIXMEMORYPOOL_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! A memory pool datatype that can be reused for matrix matrix multiplication
//! this is to prevent excessive alloc/dealloc.
class PMatrixMemoryPool {
public:
  //! Standard constructor.
  PMatrixMemoryPool(const Matrix_ps &Matrix);
  //! Standard destructor.
  ~PMatrixMemoryPool();

private:
  //! Pointer to underlying data.
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  PMatrixMemoryPool(const PMatrixMemoryPool &);
  //! Assignment operator, locked.
  PMatrixMemoryPool &operator=(const PMatrixMemoryPool &);

private:
  friend class Matrix_ps;
  friend class LoadBalancer;
};
} // namespace NTPoly
#endif
