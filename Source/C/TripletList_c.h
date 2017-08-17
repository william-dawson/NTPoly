#ifndef TRIPLETLIST_ch
#define TRIPLETLIST_ch

void ConstructTripletList_wrp(int *ih_this, const int *size);
void ResizeTripletList_wrp(int *ih_this, const int *size);
void AppendToTripletList_wrp(int *ih_this, const int *index_column,
                             const int *index_row, const double *point_value);
void SetTripletAt_wrp(int *ih_this, const int *index, const int *index_column,
                      const int *index_row, const double *point_value);
void GetTripletAt_wrp(const int *ih_this, const int *index, int *index_column,
                      const int *index_row, double *point_value);
void DestructTripletList_wrp(int *ih_this);
void SortTripletList_wrp(const int *ih_this, const int *matrix_size, 
                         int* h_sorted);

#endif
