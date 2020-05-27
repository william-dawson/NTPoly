#include "MatrixConversion.h"
#include "PSMatrix.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "MatrixConversion_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
void MatrixConversion::SnapMatrixToSparsityPattern(Matrix_ps &mata,
                                                   const Matrix_ps &matb) {
  SnapMatrixToSparsityPattern_wrp(mata.ih_this, matb.ih_this);
}
} // namespace NTPoly