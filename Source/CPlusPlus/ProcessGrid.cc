#include "ProcessGrid.h"
using namespace NTPoly;

extern "C" {
#include "ProcessGrid_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid(MPI_Comm world_comm, int process_rows,
                         int process_columns, int process_slices) {
  MPI_Fint temp_comm = MPI_Comm_c2f(world_comm);
  ConstructProcessGrid_wrp(ih_this, &temp_comm, &process_rows, &process_columns,
                           &process_slices);
}

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid(int process_rows, int process_columns,
                         int process_slices) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructProcessGrid_wrp(ih_this, &temp_comm, &process_rows, &process_columns,
                           &process_slices);
}

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid(MPI_Comm world_comm, int process_slices) {
  MPI_Fint temp_comm = MPI_Comm_c2f(world_comm);
  ConstructProcessGrid_onlyslice_wrp(ih_this, &temp_comm, &process_slices);
}

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid(int process_slices) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructProcessGrid_onlyslice_wrp(ih_this, &temp_comm, &process_slices);
}

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid() {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructProcessGrid_default_wrp(ih_this, &temp_comm);
}

//////////////////////////////////////////////////////////////////////////////
ProcessGrid::ProcessGrid(const ProcessGrid &old_grid) {
  CopyProcessGrid_wrp(old_grid.ih_this, ih_this);
}

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetMySlice() { return GetMySlice_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetMyColumn() { return GetMyColumn_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetMyRow() { return GetMyRow_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetNumSlices() { return GetNumSlices_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetNumColumns() { return GetNumColumns_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
int ProcessGrid::GetNumRows() { return GetNumRows_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
ProcessGrid::~ProcessGrid() { DestructProcessGrid_wrp(ih_this); }

////////////////////////////////////////////////////////////////////////////////
void ConstructGlobalProcessGrid(MPI_Comm world_comm, int process_rows,
                                int process_columns, int process_slices,
                                bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(world_comm);
  ConstructGlobalProcessGrid_wrp(&temp_comm, &process_rows, &process_columns,
                                 &process_slices, &be_verbose);
}

////////////////////////////////////////////////////////////////////////////////
void ConstructGlobalProcessGrid(int process_rows, int process_columns,
                                int process_slices, bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructGlobalProcessGrid_wrp(&temp_comm, &process_rows, &process_columns,
                                 &process_slices, &be_verbose);
}

////////////////////////////////////////////////////////////////////////////////
void ConstructGlobalProcessGrid(MPI_Comm world_comm, int process_slices,
                                bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(world_comm);
  ConstructGlobalProcessGrid_onlyslice_wrp(&temp_comm, &process_slices,
                                           &be_verbose);
}

////////////////////////////////////////////////////////////////////////////////
void ConstructGlobalProcessGrid(int process_slices, bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructGlobalProcessGrid_onlyslice_wrp(&temp_comm, &process_slices,
                                           &be_verbose);
}

////////////////////////////////////////////////////////////////////////////////
void ConstructGlobalProcessGrid(bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructGlobalProcessGrid_default_wrp(&temp_comm, &be_verbose);
}

////////////////////////////////////////////////////////////////////////////////
int GetGlobalMySlice() { return GetGlobalMySlice_wrp(); }

////////////////////////////////////////////////////////////////////////////////
int GetGlobalMyColumn() { return GetGlobalMyColumn_wrp(); }

////////////////////////////////////////////////////////////////////////////////
int GetGlobalMyRow() { return GetGlobalMyRow_wrp(); }

////////////////////////////////////////////////////////////////////////////////
int GetGlobalNumSlices() { return GetGlobalNumSlices_wrp(); }

////////////////////////////////////////////////////////////////////////////////
int GetGlobalNumColumns() { return GetGlobalNumColumns_wrp(); }

////////////////////////////////////////////////////////////////////////////////
int GetGlobalNumRows() { return GetGlobalNumRows_wrp(); }

////////////////////////////////////////////////////////////////////////////////
void DestructGlobalProcessGrid() { DestructGlobalProcessGrid_wrp(); }
} // namespace NTPoly
