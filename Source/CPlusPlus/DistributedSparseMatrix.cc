#include "DistributedMatrixMemoryPool.h"
#include "DistributedSparseMatrix.h"
#include "Permutation.h"
#include "TripletList.h"
using std::string;
#include <iostream>
using std::cout;
using std::endl;
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "DistributedSparseMatrix_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::DistributedSparseMatrix(int matrix_dimension) {
  ConstructEmptyMatrix_ps_wrp(ih_this, &matrix_dimension);
}

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::DistributedSparseMatrix(std::string file_name,
                                                 bool is_binary) {
  int string_length = file_name.length();
  if (is_binary) {
    ConstructMatrixFromBinary_ps_wrp(ih_this, &file_name.c_str()[0],
                                     &string_length);
  } else {
    ConstructMatrixFromMatrixMarket_ps_wrp(ih_this, &file_name.c_str()[0],
                                           &string_length);
  }
}

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::DistributedSparseMatrix(
    const DistributedSparseMatrix &matB) {
  // Constructing empty here is important because the call to Empty_wrp
  // also allocates a handle.
  int matrix_dimension = matB.GetActualDimension();
  ConstructEmptyMatrix_ps_wrp(ih_this, &matrix_dimension);
  CopyMatrix_ps_wrp(matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::WriteToBinary(std::string file_name) const {
  int string_length = file_name.length();
  WriteMatrixToBinary_ps_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  WriteMatrixToMatrixMarket_ps_wrp(ih_this, &file_name.c_str()[0],
                                   &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillFromTripletList(
    const TripletList_r &triplet_list) {
  FillMatrixFromTripletList_psr_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillFromTripletList(
    const TripletList_c &triplet_list) {
  FillMatrixFromTripletList_psc_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillDistributedPermutation(const Permutation &lb,
                                                         bool permuterows) {
  FillMatrixPermutation_ps_wrp(ih_this, lb.ih_this, &permuterows);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillIdentity() {
  FillMatrixIdentity_ps_wrp(ih_this);
}

//////////////////////////////////////////////////////////////////////////////
int DistributedSparseMatrix::GetLogicalDimension() const {
  int temp;
  GetMatrixLogicalDimension_ps_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
int DistributedSparseMatrix::GetActualDimension() const {
  int temp;
  GetMatrixActualDimension_ps_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::GetTripletList(TripletList_r &triplet_list) {
  GetMatrixTripletList_psr_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::GetTripletList(TripletList_c &triplet_list) {
  GetMatrixTripletList_psc_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::GetMatrixBlock(TripletList_r &triplet_list,
                                             int start_row, int end_row,
                                             int start_column, int end_column) {
  GetMatrixBlock_psr_wrp(ih_this, triplet_list.ih_this, &start_row, &end_row,
                        &start_column, &end_column);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::GetMatrixBlock(TripletList_c &triplet_list,
                                             int start_row, int end_row,
                                             int start_column, int end_column) {
  GetMatrixBlock_psc_wrp(ih_this, triplet_list.ih_this, &start_row, &end_row,
                        &start_column, &end_column);
}

////////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Transpose(const DistributedSparseMatrix &matA) {
  TransposeMatrix_ps_wrp(matA.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Dot(const DistributedSparseMatrix &matB) {
  return DotMatrix_ps_wrp(ih_this, matB.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Increment(const DistributedSparseMatrix &matB,
                                        double alpha, double threshold) {
  IncrementMatrix_ps_wrp(matB.ih_this, ih_this, &alpha, &threshold);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::PairwiseMultiply(
    const DistributedSparseMatrix &matA, const DistributedSparseMatrix &matB) {
  MatrixPairwiseMultiply_ps_wrp(matA.ih_this, matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Gemm(const DistributedSparseMatrix &matA,
                                   const DistributedSparseMatrix &matB,
                                   DistributedMatrixMemoryPool &memory_pool,
                                   double alpha, double beta,
                                   double threshold) {
  MatrixMultiply_ps_wrp(matA.ih_this, matB.ih_this, ih_this, &alpha, &beta,
                        &threshold, memory_pool.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Scale(double constant) {
  ScaleMatrix_ps_wrp(ih_this, &constant);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Norm() const {
  return MatrixNorm_ps_wrp(ih_this);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Trace() const {
  return MatrixTrace_ps_wrp(ih_this);
}

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::~DistributedSparseMatrix() {
  DestructMatrix_ps_wrp(ih_this);
}
} // namespace NTPoly
