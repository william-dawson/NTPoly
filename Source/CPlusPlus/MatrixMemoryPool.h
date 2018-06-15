#ifndef MATRIXMEMORYPOOL_h
#define MATRIXMEMORYPOOL_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class SparseMatrix;
////////////////////////////////////////////////////////////////////////////////
//! A memory pool datatype that can be reused for matrix matrix multiplication
//! this is to prevent excessive alloc/dealloc.
class MatrixMemoryPool_r {
public:
  //! Constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  MatrixMemoryPool_r(int columns, int rows);
  //! Destructor
  ~MatrixMemoryPool_r();

private:
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  MatrixMemoryPool_r(const MatrixMemoryPool_r &);
  //! Assignment operator, locked.
  MatrixMemoryPool_r &operator=(const MatrixMemoryPool_r &);

private:
  friend class SparseMatrix_r;
};
//! A memory pool datatype that can be reused for matrix matrix multiplication
//! this is to prevent excessive alloc/dealloc.
class MatrixMemoryPool_c {
public:
  //! Constructor.
  //!\param columns number of columns for the matrix.
  //!\param rows number of rows for the matrix.
  MatrixMemoryPool_c(int columns, int rows);
  //! Destructor
  ~MatrixMemoryPool_c();

private:
  int ih_this[SIZE_wrp];

private:
  //! Copy constructor, locked.
  MatrixMemoryPool_c(const MatrixMemoryPool_c &);
  //! Assignment operator, locked.
  MatrixMemoryPool_c &operator=(const MatrixMemoryPool_c &);

private:
  friend class SparseMatrix_c;
};
} // namespace NTPoly
#endif
