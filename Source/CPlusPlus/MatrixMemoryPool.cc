#include "MatrixMemoryPool.h"
using namespace NTPoly;

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
MatrixMemoryPool<double>::~MatrixMemoryPool() {
  DestructMatrixMemoryPool_lr_wrp(ih_this);
}
} // namespace NTPoly
