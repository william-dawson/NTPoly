#include "DistributedMatrixMemoryPool.h"
#include "DistributedSparseMatrix.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "DistributedMatrixMemoryPool_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
DistributedMatrixMemoryPool::DistributedMatrixMemoryPool(
    const DistributedSparseMatrix &Matrix) {
  ConstructDistributedMatrixMemoryPool_wrp(ih_this, Matrix.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
DistributedMatrixMemoryPool::~DistributedMatrixMemoryPool() {
  DestructDistributedMatrixMemoryPool_wrp(ih_this);
}
} // namespace NTPoly
