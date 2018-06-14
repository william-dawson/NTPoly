#ifndef MATRIXMEMORYPOOL_h
#define MATRIXMEMORYPOOL_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
template<class T> class SparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A memory pool datatype that can be reused for matrix matrix multiplication
//! this is to prevent excessive alloc/dealloc.
template <class T> class MatrixMemoryPool {
public:
  //! Constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  MatrixMemoryPool(int columns, int rows);
  //! Destructor
  ~MatrixMemoryPool();

private:
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  MatrixMemoryPool(const MatrixMemoryPool &);
  //! Assignment operator, locked.
  MatrixMemoryPool &operator=(const MatrixMemoryPool &);

private:
  template<class T2> friend class SparseMatrix;
};
} // namespace NTPoly
#endif
