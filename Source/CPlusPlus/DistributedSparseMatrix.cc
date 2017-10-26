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
  ConstructEmptyDistributedSparseMatrix_wrp(ih_this, &matrix_dimension);
}

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::DistributedSparseMatrix(std::string file_name,
                                                 bool is_binary) {
  int string_length = file_name.length();
  if (is_binary) {
    ConstructFromBinary_wrp(ih_this, &file_name.c_str()[0], &string_length);
  } else {
    ConstructFromMatrixMarket_wrp(ih_this, &file_name.c_str()[0],
                                  &string_length);
  }
}

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::DistributedSparseMatrix(
    const DistributedSparseMatrix &matB) {
  // Constructing empty here is important because the call to Empty_wrp
  // also allocates a handle.
  int matrix_dimension = matB.GetActualDimension();
  ConstructEmptyDistributedSparseMatrix_wrp(ih_this, &matrix_dimension);
  CopyDistributedSparseMatrix_wrp(matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::WriteToBinary(std::string file_name) const {
  int string_length = file_name.length();
  WriteToBinary_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::WriteToMatrixMarket(string file_name) const {
  int string_length = file_name.length();
  WriteToMatrixMarket_wrp(ih_this, &file_name.c_str()[0], &string_length);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillFromTripletList(
    const TripletList &triplet_list) {
  FillFromTripletList_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillDistributedPermutation(const Permutation &lb,
                                                         bool permuterows) {
  FillDistributedPermutation_wrp(ih_this, lb.ih_this, &permuterows);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::FillIdentity() {
  FillDistributedIdentity_wrp(ih_this);
}

//////////////////////////////////////////////////////////////////////////////
int DistributedSparseMatrix::GetLogicalDimension() const {
  int temp;
  GetLogicalDimension_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
int DistributedSparseMatrix::GetActualDimension() const {
  int temp;
  GetActualDimension_wrp(ih_this, &temp);
  return temp;
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::GetTripletList(TripletList &triplet_list) {
  GetTripletList_wrp(ih_this, triplet_list.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::RepartitionMatrix(TripletList &triplet_list,
                                                int start_row,
                                                int start_column) {
  RepartitionMatrix_wrp(ih_this, triplet_list.ih_this, &start_row,
                        &start_column);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Dot(const DistributedSparseMatrix &matB) {
  return DotDistributedSparseMatrix_wrp(ih_this, matB.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Increment(const DistributedSparseMatrix &matB,
                                        double alpha, double threshold) {
  IncrementDistributedSparseMatrix_wrp(matB.ih_this, ih_this, &alpha,
                                       &threshold);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::PairwiseMultiply(
    const DistributedSparseMatrix &matA, const DistributedSparseMatrix &matB) {
  DistributedPairwiseMultiply_wrp(matA.ih_this, matB.ih_this, ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Gemm(const DistributedSparseMatrix &matA,
                                   const DistributedSparseMatrix &matB,
                                   DistributedMatrixMemoryPool &memory_pool,
                                   double alpha, double beta,
                                   double threshold) {
  DistributedGemm_wrp(matA.ih_this, matB.ih_this, ih_this, &alpha, &beta,
                      &threshold, memory_pool.ih_this);
}

//////////////////////////////////////////////////////////////////////////////
void DistributedSparseMatrix::Scale(double constant) {
  ScaleDistributedSparseMatrix_wrp(ih_this, &constant);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Norm() const {
  return DistributedSparseNorm_wrp(ih_this);
}

//////////////////////////////////////////////////////////////////////////////
double DistributedSparseMatrix::Trace() const { return Trace_wrp(ih_this); }

//////////////////////////////////////////////////////////////////////////////
DistributedSparseMatrix::~DistributedSparseMatrix() {
  DestructDistributedSparseMatrix_wrp(ih_this);
}
} // namespace NTPoly
