#ifndef PROCESSGRID_ch
#define PROCESSGRID_ch

void ConstructGlobalProcessGrid_wrp(const int *world_comm,
                                    const int *process_rows,
                                    const int *process_columns,
                                    const int *process_slices,
                                    const bool *be_verbose);
void CopyProcessGrid_wrp(const int *ih_old_grid, int* ih_new_grid);
int GetGlobalMySlice_wrp();
int GetGlobalMyColumn_wrp();
int GetGlobalMyRow_wrp();
void DestructGlobalProcessGrid_wrp();

void ConstructProcessGrid_wrp(int *ih_grid, const int *world_comm,
                              const int *process_rows,
                              const int *process_columns,
                              const int *process_slices,
                              const bool *be_verbose);
int GetMySlice_wrp(const int *ih_grid);
int GetMyColumn_wrp(const int *ih_grid);
int GetMyRow_wrp(const int *ih_grid);
void DestructProcessGrid_wrp(int *ih_grid);

#endif
