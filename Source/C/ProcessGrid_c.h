#ifndef PROCESSGRID_ch
#define PROCESSGRID_ch

void ConstructProcessGrid_wrp(const int *world_comm, const int *process_rows,
                              const int *process_columns,
                              const int *process_slices,
                              const bool *be_verbose);
int GetMySlice_wrp();
int GetMyColumn_wrp();
int GetMyRow_wrp();
void DestructProcessGrid_wrp();

#endif
