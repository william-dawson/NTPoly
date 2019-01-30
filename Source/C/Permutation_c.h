#ifndef PERMUTATION_ch
#define PERMUTATION_ch

void ConstructDefaultPermutation_wrp(int *ih_this, const int *matrix_dimension);
void ConstructReversePermutation_wrp(int *ih_this, const int *matrix_dimension);
void ConstructRandomPermutation_wrp(int *ih_this, const int *matrix_dimension,
                                    const int *seed);
void DestructPermutation_wrp(int *ih_this);

#endif
