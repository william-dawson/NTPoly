#ifndef TRIPLETLIST_ch
#define TRIPLETLIST_ch

void ConstructTripletList_r_wrp(int *ih_this, const int *size);
void ResizeTripletList_r_wrp(int *ih_this, const int *size);
void AppendToTripletList_r_wrp(int *ih_this, const int *index_column,
                               const int *index_row, const double *point_value);
void SetTripletAt_r_wrp(int *ih_this, const int *index, const int *index_column,
                        const int *index_row, const double *point_value);
void GetTripletAt_r_wrp(const int *ih_this, const int *index, int *index_column,
                        const int *index_row, double *point_value);
void DestructTripletList_r_wrp(int *ih_this);
void SortTripletList_r_wrp(const int *ih_this, const int *matrix_size,
                           int *h_sorted);
int GetTripletListSize_r_wrp(const int *ih_this);

void ConstructTripletList_c_wrp(int *ih_this, const int *size);
void ResizeTripletList_c_wrp(int *ih_this, const int *size);
void AppendToTripletList_c_wrp(int *ih_this, const int *index_column,
                               const int *index_row, const double *point_value);
void SetTripletAt_c_wrp(int *ih_this, const int *index, const int *index_column,
                        const int *index_row, const double *point_value);
void GetTripletAt_c_wrp(const int *ih_this, const int *index, int *index_column,
                        const int *index_row, double *point_value);
void DestructTripletList_c_wrp(int *ih_this);
void SortTripletList_c_wrp(const int *ih_this, const int *matrix_size,
                           int *h_sorted);
int GetTripletListSize_c_wrp(const int *ih_this);

#endif
