#ifndef DISTRIBUTEDSPARSEMATRIX_ch
#define DISTRIBUTEDSPARSEMATRIX_ch

void ConstructEmptyMatrix_ps_wrp(int *ih_this, const int *matrix_dim);
void ConstructEmptyMatrixPG_ps_wrp(int *ih_this, const int *matrix_dim,
                                   const int *ih_grid);
void CopyMatrix_ps_wrp(const int *ih_matA, int *ih_matB);
void DestructMatrix_ps_wrp(int *ih_this);
void ConstructMatrixFromMatrixMarket_ps_wrp(int *ih_this, const char *file_name,
                                            const int *name_size);
void ConstructMatrixFromBinary_ps_wrp(int *ih_this, const char *file_name,
                                      const int *name_size);
void ConstructMatrixFromMatrixMarketPG_ps_wrp(int *ih_this,
                                              const char *file_name,
                                              const int *name_size,
                                              const int *ih_grid);
void ConstructMatrixFromBinaryPG_ps_wrp(int *ih_this, const char *file_name,
                                        const int *name_size,
                                        const int *ih_grid);
void WriteMatrixToBinary_ps_wrp(const int *ih_this, const char *file_name,
                                const int *name_size);
void WriteMatrixToMatrixMarket_ps_wrp(const int *ih_this, const char *file_name,
                                      const int *name_size);
void FillMatrixFromTripletList_psr_wrp(const int *ih_this,
                                       const int *ih_triplet_list);
void FillMatrixFromTripletList_psc_wrp(const int *ih_this,
                                       const int *ih_triplet_list);
void FillMatrixPermutation_ps_wrp(int *ih_this, const int *ih_permutation,
                                  const bool *permuterows);
void FillMatrixIdentity_ps_wrp(int *ih_this);
void GetMatrixActualDimension_ps_wrp(const int *ih_this, int *size);
void GetMatrixLogicalDimension_ps_wrp(const int *ih_this, int *size);
void GetMatrixSize_ps_wrp(const int *ih_this, long int *size);
void GetMatrixTripletList_psr_wrp(const int *ih_this, int *ih_triplet_list);
void GetMatrixTripletList_psc_wrp(const int *ih_this, int *ih_triplet_list);
void GetMatrixBlock_psr_wrp(const int *ih_this, int *ih_triplet_list,
                            int *start_row, int *end_row, int *start_column,
                            int *end_column);
void GetMatrixBlock_psc_wrp(const int *ih_this, int *ih_triplet_list,
                            int *start_row, int *end_row, int *start_column,
                            int *end_column);
void GetMatrixSlice_wrp(const int *ih_this, int *ih_submatrix, int *start_row,
                        int *end_row, int *start_column, int *end_column);
void TransposeMatrix_ps_wrp(const int *ih_matA, int *ih_transmat);
void ConjugateMatrix_ps_wrp(int *ih_matA);
void ConjugateMatrix_ps_wrp(int *ih_matA);
void ResizeMatrix_ps_wrp(int *ih_this, const int *new_size);
void GetMatrixProcessGrid_ps_wrp(const int *ih_this, int *ih_grid);

void DotMatrix_psr_wrp(const int *ih_matA, const int *ih_matB, double *product);
void DotMatrix_psc_wrp(const int *ih_matA, const int *ih_matB,
                       double *product_real, double *product_imag);
void IncrementMatrix_ps_wrp(const int *ih_matA, int *ih_matB,
                            const double *alpha_in, const double *threshold_in);
void MatrixPairwiseMultiply_ps_wrp(const int *ih_matA, const int *ih_matB,
                                   int *ih_matC);
void MatrixMultiply_ps_wrp(const int *ih_matA, const int *ih_matB, int *ih_matC,
                           const double *alpha_in, const double *beta_in,
                           const double *threshold_in, int *ih_memory_pool_in);
void ScaleMatrix_ps_wrp(int *ih_this, const double *constant);
double MatrixNorm_ps_wrp(const int *ih_this);
void MatrixTrace_ps_wrp(const int *ih_this, double *trace_val);

#endif
