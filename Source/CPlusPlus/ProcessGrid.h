#ifndef PROCESSGRID_h
#define PROCESSGRID_h
#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
//! Construct the process grid.
//! \param[in] world_comm a communicator that every process in the grid is
//! a part of.
//! \param[in] process_rows number of grid rows.
//! \param[in] process_columns number of grid columns.
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructProcessGrid(MPI_Comm world_comm, int process_rows,
                          int process_columns, int process_slices,
                          bool be_verbose = false);
//! Construct the process grid from comm world
//! \param[in] process_rows number of grid rows.
//! \param[in] process_columns number of grid columns.
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructProcessGrid(int process_rows, int process_columns,
                          int process_slices, bool be_verbose = false);
//! Get the slice of the current process.
int GetMySlice();
//! Get the column of the current process.
int GetMyColumn();
//! Get the row of the current process.
int GetMyRow();
} // namespace NTPoly
#endif
