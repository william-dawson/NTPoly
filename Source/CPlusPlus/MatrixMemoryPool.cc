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

MatrixMemoryPool_r::MatrixMemoryPool_r(int columns, int rows) {
  ConstructMatrixMemoryPool_lr_wrp(ih_this, &rows, &columns);
}

MatrixMemoryPool_c::MatrixMemoryPool_c(int columns, int rows) {
  ConstructMatrixMemoryPool_lc_wrp(ih_this, &rows, &columns);
}

////////////////////////////////////////////////////////////////////////////////

MatrixMemoryPool_r::~MatrixMemoryPool_r() {
  DestructMatrixMemoryPool_lc_wrp(ih_this);
}

MatrixMemoryPool_c::~MatrixMemoryPool_c() {
  DestructMatrixMemoryPool_lc_wrp(ih_this);
}
} // namespace NTPoly
