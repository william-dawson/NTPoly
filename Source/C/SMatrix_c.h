#ifndef SparseMatrix_ch
#define SparseMatrix_ch

void ConstructMatrixFromFile_lsr_wrp(int *ih_this, const char *file_name,
                                     const int *name_size);
void ConstructMatrixFromTripletList_lsr_wrp(int *ih_this,
                                            const int *ih_triplet_list,
                                            const int *rows,
                                            const int *columns);
void ConstructZeroMatrix_lsr_wrp(int *ih_this, const int *rows,
                                 const int *columns);
void DestructMatrix_lsr_wrp(int *ih_this);
void CopyMatrix_lsr_wrp(const int *ih_matA, int *ih_matB);
void GetMatrixRows_lsr_wrp(const int *ih_this, int *rows);
void GetMatrixColumns_lsr_wrp(const int *ih_this, int *columns);
void ExtractMatrixRow_lsr_wrp(const int *ih_this, int *row_number,
                              int *ih_row_out);
void ExtractMatrixColumn_lsr_wrp(const int *ih_this, int *column_number,
                                 int *ih_column_out);
void ScaleMatrix_lsr_wrp(int *ih_this, const double *constant);
void IncrementMatrix_lsr_wrp(const int *ih_matA, int *ih_matB,
                             const double *alpha_in,
                             const double *threshold_in);
void DotMatrix_lsr_wrp(const int *ih_matA, const int *ih_matB, double *product);
void PairwiseMultiplyMatrix_lsr_wrp(const int *ih_matA, const int *ih_matB,
                                    int *ih_matC);
void MatrixMultiply_lsr_wrp(const int *ih_matA, const int *ih_matB,
                            int *ih_matC, const bool *IsATransposed,
                            const bool *IsBTransposed, const double *alpha,
                            const double *beta, const double *threshold,
                            int *ih_matrix_memory_pool);
void TransposeMatrix_lsr_wrp(const int *ih_matA, int *ih_matAT);
void PrintMatrix_lsr_wrp(const int *ih_this);
void PrintMatrixF_lsr_wrp(const int *ih_this, const char *file_name,
                          const int *name_size);
void MatrixToTripletList_lsr_wrp(const int *ih_this, int *ih_triplet_list);

void ConstructMatrixFromFile_lsc_wrp(int *ih_this, const char *file_name,
                                     const int *name_size);
void ConstructMatrixFromTripletList_lsc_wrp(int *ih_this,
                                            const int *ih_triplet_list,
                                            const int *rows,
                                            const int *columns);
void ConstructZeroMatrix_lsc_wrp(int *ih_this, const int *rows,
                                 const int *columns);
void DestructMatrix_lsc_wrp(int *ih_this);
void CopyMatrix_lsc_wrp(const int *ih_matA, int *ih_matB);
void GetMatrixRows_lsc_wrp(const int *ih_this, int *rows);
void GetMatrixColumns_lsc_wrp(const int *ih_this, int *columns);
void ExtractMatrixRow_lsc_wrp(const int *ih_this, int *row_number,
                              int *ih_row_out);
void ExtractMatrixColumn_lsc_wrp(const int *ih_this, int *column_number,
                                 int *ih_column_out);
void ScaleMatrix_lsc_wrp(int *ih_this, const double *constant);
void IncrementMatrix_lsc_wrp(const int *ih_matA, int *ih_matB,
                             const double *alpha_in,
                             const double *threshold_in);
void DotMatrix_lsc_wrp(const int *ih_matA, const int *ih_matB,
                       double *product_real, double *product_complex);
void PairwiseMultiplyMatrix_lsc_wrp(const int *ih_matA, const int *ih_matB,
                                    int *ih_matC);
void MatrixMultiply_lsc_wrp(const int *ih_matA, const int *ih_matB,
                            int *ih_matC, const bool *IsATransposed,
                            const bool *IsBTransposed, const double *alpha,
                            const double *beta, const double *threshold,
                            int *ih_matrix_memory_pool);
void TransposeMatrix_lsc_wrp(const int *ih_matA, int *ih_matAT);
void ConjugateMatrix_lsc_wrp(int *ih_matA);
void PrintMatrix_lsc_wrp(const int *ih_this);
void PrintMatrixF_lsc_wrp(const int *ih_this, const char *file_name,
                          const int *name_size);
void MatrixToTripletList_lsc_wrp(const int *ih_this, int *ih_triplet_list);

#endif
