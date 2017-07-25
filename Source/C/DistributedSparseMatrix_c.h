#ifndef DISTRIBUTEDSPARSEMATRIX_ch
#define DISTRIBUTEDSPARSEMATRIX_ch

void ConstructEmpty_wrp(int *ih_this, const int *matrix_dim);
void CopyDistributedSparseMatrix_wrp(const int *ih_matA, int *ih_matB);
void DestructDistributedSparseMatrix_wrp(int *ih_this);
void ConstructFromMatrixMarket_wrp(int *ih_this, const char *file_name,
                                   const int *name_size);
void ConstructFromBinary_wrp(int *ih_this, const char *file_name,
                             const int *name_size);
void WriteToBinary_wrp(const int *ih_this, const char *file_name,
                       const int *name_size);
void WriteToMatrixMarket_wrp(const int *ih_this, const char *file_name,
                             const int *name_size);
void FillFromTripletList_wrp(const int *ih_this, const int *ih_triplet_list);
void FillDistributedPermutation_wrp(int *ih_this, const int *ih_permutation,
                                    const bool *permuterows);
void FillDistributedIdentity_wrp(int *ih_this);
void GetActualDimension_wrp(const int *ih_this, int *size);
void GetLogicalDimension_wrp(const int *ih_this, int *size);
double DotDistributedSparseMatrix_wrp(const int *ih_matA, const int *ih_matB);
void IncrementDistributedSparseMatrix_wrp(const int *ih_matA, int *ih_matB,
                                          const double *alpha_in,
                                          const double *threshold_in);
void DistributedPairwiseMultiply_wrp(const int *ih_matA, const int *ih_matB,
                                     int *ih_matC);
void DistributedGemm_wrp(const int *ih_matA, const int *ih_matB, int *ih_matC,
                         const double *alpha_in, const double *beta_in,
                         const double *threshold_in, int *ih_memory_pool_in);
void ScaleDistributedSparseMatrix_wrp(int *ih_this, const double *constant);
double DistributedSparseNorm_wrp(const int *ih_this);
double Trace_wrp(const int *ih_this);

#endif
