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
SparseMatrix_r::SparseMatrix_r(int columns, int rows) {
  ConstructZeroMatrix_lsr_wrp(ih_this, &rows, &columns);
}
SparseMatrix_c::SparseMatrix_c(int columns, int rows) {
  ConstructZeroMatrix_lsc_wrp(ih_this, &rows, &columns);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix_r::SparseMatrix_r(std::string file_name) {
  int string_length = file_name.length();
  ConstructMatrixFromFile_lsr_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
}
SparseMatrix_c::SparseMatrix_c(std::string file_name) {
  int string_length = file_name.length();
  ConstructMatrixFromFile_lsc_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
}

////////////////////////////////////////////////////////////////////////////////

SparseMatrix_r::SparseMatrix_r(const TripletList_r &list, int rows,
                               int columns) {
  ConstructMatrixFromTripletList_lsr_wrp(ih_this, list.ih_this, &rows,
                                         &columns);
}

SparseMatrix_c::SparseMatrix_c(const TripletList_c &list, int rows,
                               int columns) {
  ConstructMatrixFromTripletList_lsc_wrp(ih_this, list.ih_this, &rows,
                                         &columns);
}

////////////////////////////////////////////////////////////////////////////////

SparseMatrix_r::SparseMatrix_r(const SparseMatrix_r &matB) {
  CopyMatrix_lsr_wrp(matB.ih_this, ih_this);
}

SparseMatrix_c::SparseMatrix_c(const SparseMatrix_c &matB) {
  CopyMatrix_lsc_wrp(matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
SparseMatrix_r::~SparseMatrix_r() { DestructMatrix_lsr_wrp(ih_this); }
SparseMatrix_c::~SparseMatrix_c() { DestructMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int SparseMatrix_r::GetRows() const {
  int temp;
  GetMatrixRows_lsr_wrp(ih_this, &temp);
  return temp;
}
int SparseMatrix_c::GetRows() const {
  int temp;
  GetMatrixRows_lsc_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
int SparseMatrix_r::GetColumns() const {
  int temp;
  GetMatrixColumns_lsr_wrp(ih_this, &temp);
  return temp;
}
int SparseMatrix_c::GetColumns() const {
  int temp;
  GetMatrixColumns_lsc_wrp(ih_this, &temp);
  return temp;
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::ExtractRow(int row_number, SparseMatrix_r &row_out) const {
  int temp = row_number + 1;
  ExtractMatrixRow_lsr_wrp(ih_this, &temp, row_out.ih_this);
}

void SparseMatrix_c::ExtractRow(int row_number, SparseMatrix_c &row_out) const {
  int temp = row_number + 1;
  ExtractMatrixRow_lsc_wrp(ih_this, &temp, row_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::ExtractColumn(int column_number,
                                   SparseMatrix_r &column_out) const {
  int temp = column_number + 1;
  ExtractMatrixColumn_lsr_wrp(ih_this, &temp, column_out.ih_this);
}

void SparseMatrix_c::ExtractColumn(int column_number,
                                   SparseMatrix_c &column_out) const {
  int temp = column_number + 1;
  ExtractMatrixColumn_lsc_wrp(ih_this, &temp, column_out.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix_r::Scale(double constant) {
  ScaleMatrix_lsr_wrp(ih_this, &constant);
}
void SparseMatrix_c::Scale(double constant) {
  ScaleMatrix_lsc_wrp(ih_this, &constant);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::Increment(const SparseMatrix_r &matB, double alpha,
                               double threshold) {
  IncrementMatrix_lsr_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

void SparseMatrix_c::Increment(const SparseMatrix_c &matB, double alpha,
                               double threshold) {
  IncrementMatrix_lsc_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

////////////////////////////////////////////////////////////////////////////////

double SparseMatrix_r::Dot(const SparseMatrix_r &matB) const {
  return DotMatrix_lsr_wrp(ih_this, matB.ih_this);
}

double SparseMatrix_c::Dot(const SparseMatrix_c &matB) const {
  return DotMatrix_lsc_wrp(ih_this, matB.ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::PairwiseMultiply(const NTPoly::SparseMatrix_r &matA,
                                      const NTPoly::SparseMatrix_r &matB) {
  PairwiseMultiplyMatrix_lsr_wrp(matA.ih_this, matB.ih_this, ih_this);
}

void SparseMatrix_c::PairwiseMultiply(const NTPoly::SparseMatrix_c &matA,
                                      const NTPoly::SparseMatrix_c &matB) {
  PairwiseMultiplyMatrix_lsc_wrp(matA.ih_this, matB.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::Gemm(const SparseMatrix_r &matA,
                          const SparseMatrix_r &matB, bool isATransposed,
                          bool isBTransposed, double alpha, double beta,
                          double threshold, MatrixMemoryPool_r &memory_pool) {
  MatrixMultiply_lsr_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed,
                         &isBTransposed, &alpha, &beta, &threshold,
                         memory_pool.ih_this);
}

void SparseMatrix_c::Gemm(const SparseMatrix_c &matA,
                          const SparseMatrix_c &matB, bool isATransposed,
                          bool isBTransposed, double alpha, double beta,
                          double threshold, MatrixMemoryPool_c &memory_pool) {
  MatrixMultiply_lsc_wrp(matA.ih_this, matB.ih_this, ih_this, &isATransposed,
                         &isBTransposed, &alpha, &beta, &threshold,
                         memory_pool.ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix_r::Transpose(const SparseMatrix_r &matA) {
  TransposeMatrix_lsr_wrp(matA.ih_this, ih_this);
}

void SparseMatrix_c::Transpose(const SparseMatrix_c &matA) {
  TransposeMatrix_lsc_wrp(matA.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix_c::Conjugate() { ConjugateMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix_r::Print() const { PrintMatrix_lsr_wrp(ih_this); }
void SparseMatrix_c::Print() const { PrintMatrix_lsc_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void SparseMatrix_r::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  PrintMatrixF_lsr_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

void SparseMatrix_c::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  PrintMatrixF_lsc_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::MatrixToTripletList(TripletList_r &triplet_list) const {
  MatrixToTripletList_lsr_wrp(ih_this, triplet_list.ih_this);
}

void SparseMatrix_c::MatrixToTripletList(TripletList_c &triplet_list) const {
  MatrixToTripletList_lsc_wrp(ih_this, triplet_list.ih_this);
}

////////////////////////////////////////////////////////////////////////////////

void SparseMatrix_r::EigenDecomposition(NTPoly::SparseMatrix_r &MatV,
                                        double threshold) const {
  EigenDecomposition_lsr_wrp(ih_this, MatV.ih_this, &threshold);
}

void SparseMatrix_c::EigenDecomposition(NTPoly::SparseMatrix_c &MatV,
                                        double threshold) const {
  EigenDecomposition_lsc_wrp(ih_this, MatV.ih_this, &threshold);
}
} // namespace NTPoly
