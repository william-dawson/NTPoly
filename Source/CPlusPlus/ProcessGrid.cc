#include "ProcessGrid.h"
using namespace NTPoly;

extern "C" {
#include "ProcessGrid_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
void ConstructProcessGrid(MPI_Comm world_comm, int process_rows,
                          int process_columns, int process_slices,
                          bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(world_comm);
  ConstructProcessGrid_wrp(&temp_comm, &process_rows, &process_columns,
                           &process_slices, &be_verbose);
}
void ConstructProcessGrid(int process_rows, int process_columns,
                          int process_slices, bool be_verbose) {
  MPI_Fint temp_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  ConstructProcessGrid_wrp(&temp_comm, &process_rows, &process_columns,
                           &process_slices, &be_verbose);
}
int GetMySlice() { return GetMySlice_wrp(); }
int GetMyColumn() { return GetMyColumn_wrp(); }
int GetMyRow() { return GetMyRow_wrp(); }
} // namespace NTPoly
