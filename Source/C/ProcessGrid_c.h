#ifndef PROCESSGRID_ch
#define PROCESSGRID_ch

void ConstructGlobalProcessGrid_wrp(const int *world_comm,
                                    const int *process_rows,
                                    const int *process_columns,
                                    const int *process_slices,
                                    const bool *be_verbose);
void ConstructGlobalProcessGrid_onlyslice_wrp(const int *world_comm,
                                              const int *process_slices,
                                              const bool *be_verbose);
void ConstructGlobalProcessGrid_default_wrp(const int *world_comm,
                                            const bool *be_verbose);
void CopyProcessGrid_wrp(const int *ih_old_grid, int *ih_new_grid);
int GetGlobalMySlice_wrp();
int GetGlobalMyColumn_wrp();
int GetGlobalMyRow_wrp();
bool GetGlobalIsRoot_wrp();
int GetGlobalNumSlices_wrp();
int GetGlobalNumColumns_wrp();
int GetGlobalNumRows_wrp();
void WriteGlobalProcessGridInfo_wrp();
void DestructGlobalProcessGrid_wrp();

void ConstructProcessGrid_wrp(int *ih_grid, const int *world_comm,
                              const int *process_rows,
                              const int *process_columns,
                              const int *process_slices);
void ConstructProcessGrid_onlyslice_wrp(int *ih_grid, const int *world_comm,
                                        const int *process_slices);
void ConstructProcessGrid_default_wrp(int *ih_grid, const int *world_comm);
int GetMySlice_wrp(const int *ih_grid);
int GetMyColumn_wrp(const int *ih_grid);
int GetMyRow_wrp(const int *ih_grid);
int GetNumSlices_wrp(const int *ih_grid);
int GetNumColumns_wrp(const int *ih_grid);
int GetNumRows_wrp(const int *ih_grid);
void WriteProcessGridInfo_wrp(const int *ih_grid);
void DestructProcessGrid_wrp(int *ih_grid);

#endif
