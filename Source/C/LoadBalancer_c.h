#ifndef LOADBALANCER_ch
#define LOADBALANCER_ch

void PermuteMatrix_wrp(const int *ih_mat_in, int *ih_mat_out,
                       const int *ih_permutation, int *ih_memorypool);
void UndoPermuteMatrix_wrp(const int *ih_mat_in, int *ih_mat_out,
                           const int *ih_permutation, int *ih_memorypool);

#endif
