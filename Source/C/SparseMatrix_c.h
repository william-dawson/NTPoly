#ifndef SparseMatrix_ch
#define SparseMatrix_ch

void ConstructEmptySparseMatrix_wrp(int *ih_this, const int *columns,
                                    const int *rows);
void ConstructSparseMatrixFromFile_wrp(int *ih_this, const char *file_name,
                                       const int *name_size);
void ConstructFromTripletList_wrp(int *ih_this, const int *ih_triplet_list,
                                  const int *rows, const int *columns);
void CopySparseMatrix_wrp(const int *ih_matA, int *ih_matB);
void GetRows_wrp(const int *ih_this, int *rows);
void GetColumns_wrp(const int *ih_this, int *columns);
void ScaleSparseMatrix_wrp(int *ih_this, const double *constant);
void IncrementSparseMatrix_wrp(const int *ih_matA, int *ih_matB,
                               const double *alpha_in,
                               const double *threshold_in);
double DotSparseMatrix_wrp(const int *ih_matA, const int *ih_matB);
void PairwiseMultiplySparseMatrix_wrp(const int *ih_matA, const int *ih_matB,
                                      int *ih_matC);
void Gemm_wrp(const int *ih_matA, const int *ih_matB, int *ih_matC,
              const bool *IsATransposed, const bool *IsBTransposed,
              const double *alpha, const double *beta, const double *threshold,
              int *ih_matrix_memory_pool);
void TransposeSparseMatrix_wrp(const int *ih_matA, int *ih_matAT);
void PrintSparseMatrix_wrp(const int *ih_this);
void PrintSparseMatrixF_wrp(const int *ih_this, const char *file_name,
                            const int *name_size);
void MatrixToTripletList_wrp(const int *ih_this, int *ih_triplet_list);
void DestructSparseMatrix_wrp(int *ih_this);

#endif
