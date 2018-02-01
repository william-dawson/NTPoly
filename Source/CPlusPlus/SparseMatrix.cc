#include "MatrixMemoryPool.h"
#include "SparseMatrix.h"
#include "TripletList.h"
using std::string;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "SparseMatrix_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(int columns, int rows) {
  ConstructZeroSparseMatrix_wrp(ih_this, &columns, &rows);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(std::string file_name) {
  int string_length = file_name.length();
  ConstructSparseMatrixFromFile_wrp(ih_this, &file_name.c_str()[0],
                                    &string_length);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(const TripletList &list, int rows, int columns) {
  ConstructFromTripletList_wrp(ih_this, list.ih_this, &rows, &columns);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix::SparseMatrix(const SparseMatrix &matB) {
  CopySparseMatrix_wrp(matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix::~SparseMatrix() { DestructSparseMatrix_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int SparseMatrix::GetRows() const {
  int temp;
  GetRows_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
int SparseMatrix::GetColumns() const {
  int temp;
  GetColumns_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::ExtractRow(int row_number, SparseMatrix &row_out) const {
  int temp = row_number + 1;
  ExtractRow_wrp(ih_this, &temp, row_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::ExtractColumn(int column_number,
                                 SparseMatrix &column_out) const {
  int temp = column_number + 1;
  ExtractColumn_wrp(ih_this, &temp, column_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Scale(double constant) {
  ScaleSparseMatrix_wrp(ih_this, &constant);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Increment(const SparseMatrix &matB, double alpha,
                             double threshold) {
  IncrementSparseMatrix_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

////////////////////////////////////////////////////////////////////////////////
double SparseMatrix::Dot(const SparseMatrix &matB) const {
  return DotSparseMatrix_wrp(ih_this, matB.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void SparseMatrix::PairwiseMultiply(const NTPoly::SparseMatrix &matA,
                                    const NTPoly::SparseMatrix &matB) {
  PairwiseMultiplySparseMatrix_wrp(matA.ih_this, matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Gemm(const SparseMatrix &matA, const SparseMatrix &matB,
                        bool isATransposed, bool isBTransposed, double alpha,
                        double beta, double threshold,
                        MatrixMemoryPool &memory_pool) {
  Gemm_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed, &isBTransposed,
           &alpha, &beta, &threshold, memory_pool.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Transpose(const SparseMatrix &matA) {
  TransposeSparseMatrix_wrp(matA.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::Print() { PrintSparseMatrix_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::WriteToMatrixMarket(string file_name) {
  int string_length = file_name.length();
  PrintSparseMatrixF_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix::MatrixToTripletList(TripletList &triplet_list) {
  MatrixToTripletList_wrp(ih_this, triplet_list.ih_this);
}
} // namespace NTPoly
