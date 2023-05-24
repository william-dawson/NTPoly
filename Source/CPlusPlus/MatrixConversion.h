#ifndef MATRIXCONVERSION_h
#define MATRIXCONVERSION_h

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
////////////////////////////////////////////////////////////////////////////////
//! Helper routines for converting an NTPoly matrix to other program types.
class MatrixConversion {
public:
  //! Some codes use a fixed sparsity pattern for a matrix instead of filtering
  //! small values. Using this routine, the matrix is filled to have the same
  //! pattern as the second matrix argument. Zeros of the sparsity pattern are
  //! left in, whereas values outside the sparsity are removed. This can
  //! faciliate conversion between formats.
  //! \param mata the matrix to modify.
  //! \param matb the matrix which defines the sparsity pattern.
  static void SnapMatrixToSparsityPattern(Matrix_ps &mata,
                                          const Matrix_ps &matb);
};
} // namespace NTPoly
#endif
