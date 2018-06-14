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
template <> SparseMatrix<double>::SparseMatrix(int columns, int rows) {
  ConstructZeroMatrix_lsr_wrp(ih_this, &columns, &rows);
}

////////////////////////////////////////////////////////////////////////////////
template <> SparseMatrix<double>::SparseMatrix(std::string file_name) {
  int string_length = file_name.length();
  ConstructMatrixFromFile_lsr_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
}

////////////////////////////////////////////////////////////////////////////////
template <>
SparseMatrix<double>::SparseMatrix(const TripletList<double> &list, int rows,
                                   int columns) {
  ConstructMatrixFromTripletList_lsr_wrp(ih_this, list.ih_this, &rows,
                                         &columns);
}

////////////////////////////////////////////////////////////////////////////////
template <>
SparseMatrix<double>::SparseMatrix(const SparseMatrix<double> &matB) {
  CopyMatrix_lsr_wrp(matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <> SparseMatrix<double>::~SparseMatrix() {
  DestructMatrix_lsr_wrp(ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <> int SparseMatrix<double>::GetRows() const {
  int temp;
  GetMatrixRows_lsr_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
template <> int SparseMatrix<double>::GetColumns() const {
  int temp;
  GetMatrixColumns_lsr_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::ExtractRow(int row_number,
                                      SparseMatrix<double> &row_out) const {
  int temp = row_number + 1;
  ExtractMatrixRow_lsr_wrp(ih_this, &temp, row_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::ExtractColumn(
    int column_number, SparseMatrix<double> &column_out) const {
  int temp = column_number + 1;
  ExtractMatrixColumn_lsr_wrp(ih_this, &temp, column_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <> void SparseMatrix<double>::Scale(double constant) {
  ScaleMatrix_lsr_wrp(ih_this, &constant);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::Increment(const SparseMatrix<double> &matB,
                                     double alpha, double threshold) {
  IncrementMatrix_lsr_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

////////////////////////////////////////////////////////////////////////////////
template <>
double SparseMatrix<double>::Dot(const SparseMatrix<double> &matB) const {
  return DotMatrix_lsr_wrp(ih_this, matB.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::PairwiseMultiply(
    const NTPoly::SparseMatrix<double> &matA,
    const NTPoly::SparseMatrix<double> &matB) {
  PairwiseMultiplyMatrix_lsr_wrp(matA.ih_this, matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::Gemm(const SparseMatrix<double> &matA,
                                const SparseMatrix<double> &matB,
                                bool isATransposed, bool isBTransposed,
                                double alpha, double beta, double threshold,
                                MatrixMemoryPool<double> &memory_pool) {
  MatrixMultiply_lsr_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed,
                         &isBTransposed, &alpha, &beta, &threshold,
                         memory_pool.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::Transpose(const SparseMatrix<double> &matA) {
  TransposeMatrix_lsr_wrp(matA.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <> void SparseMatrix<double>::Print() { PrintMatrix_lsr_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
template <> void SparseMatrix<double>::WriteToMatrixMarket(string file_name) {
  int string_length = file_name.length();
  PrintMatrixF_lsr_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::MatrixToTripletList(
    TripletList<double> &triplet_list) {
  MatrixToTripletList_lsr_wrp(ih_this, triplet_list.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
template <>
void SparseMatrix<double>::EigenDecomposition(
    NTPoly::SparseMatrix<double> &MatV, double threshold) {
  EigenDecomposition_lsr_wrp(ih_this, MatV.ih_this, &threshold);
}
} // namespace NTPoly
