#include "MatrixMemoryPool.h"
using namespace NTPoly;

#include <complex.h>

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "MatrixMemoryPool_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
template<>
MatrixMemoryPool<double>::MatrixMemoryPool(int columns, int rows) {
  ConstructMatrixMemoryPool_lr_wrp(ih_this, &columns, &rows);
}
template<>
MatrixMemoryPool<double _Complex>::MatrixMemoryPool(int columns, int rows) {
  ConstructMatrixMemoryPool_lc_wrp(ih_this, &columns, &rows);
}

////////////////////////////////////////////////////////////////////////////////
template<>
MatrixMemoryPool<double>::~MatrixMemoryPool() {
  DestructMatrixMemoryPool_lc_wrp(ih_this);
}
template<>
MatrixMemoryPool<double _Complex>::~MatrixMemoryPool() {
  DestructMatrixMemoryPool_lc_wrp(ih_this);
}
} // namespace NTPoly
