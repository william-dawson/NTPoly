#ifndef PROCESSGRID_h
#define PROCESSGRID_h
#include <mpi.h>

#include "Wrapper.h"

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
class Matrix_ps;
class ProcessGrid {
public:
  //! Construct the process grid.
  //! \param[in] world_comm a communicator that every process in the grid is
  //! a part of.
  //! \param[in] process_rows number of grid rows.
  //! \param[in] process_columns number of grid columns.
  //! \param[in] process_slices number of grid slices.
  ProcessGrid(MPI_Comm world_comm, int process_rows, int process_columns,
              int process_slices);
  //! Construct the process grid from comm world
  //! \param[in] process_rows number of grid rows.
  //! \param[in] process_columns number of grid columns.
  //! \param[in] process_slices number of grid slices.
  ProcessGrid(int process_rows, int process_columns, int process_slices);
  //! Construct the process grid.
  //! \param[in] world_comm a communicator that every process in the grid is
  //! a part of.
  //! \param[in] process_rows number of grid rows.
  ProcessGrid(MPI_Comm world_comm, int process_slices);
  //! Construct the process grid from comm world
  //! \param[in] process_slices number of grid slices.
  ProcessGrid(int process_slices);
  //! Copy constructor.
  //!\param old_grid to copy from.
  ProcessGrid(const ProcessGrid &old_grid);

public:
  //! Get the slice of the current process.
  int GetMySlice();
  //! Get the column of the current process.
  int GetMyColumn();
  //! Get the row of the current process.
  int GetMyRow();
  //! Get the number of slices in this grid.
  int GetNumSlices();
  //! Get the number of columns in this grid.
  int GetNumColumns();
  //! Get the number of rows in this grid.
  int GetNumRows();

public:
  //! Standard destructor
  ~ProcessGrid();

private:
  int ih_this[SIZE_wrp];
  //! Assignment operator, locked.
  ProcessGrid &operator=(const ProcessGrid &);
  friend class Matrix_ps;
};
////////////////////////////////////////////////////////////////////////////////
//! Construct the global process grid.
//! \param[in] world_comm a communicator that every process in the grid is
//! a part of.
//! \param[in] process_rows number of grid rows.
//! \param[in] process_columns number of grid columns.
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructGlobalProcessGrid(MPI_Comm world_comm, int process_rows,
                                int process_columns, int process_slices,
                                bool be_verbose = false);
//! Construct the global process grid from comm world
//! \param[in] process_rows number of grid rows.
//! \param[in] process_columns number of grid columns.
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructGlobalProcessGrid(int process_rows, int process_columns,
                                int process_slices, bool be_verbose = false);
//! Construct the global process grid.
//! \param[in] world_comm a communicator that every process in the grid is
//! a part of.
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructGlobalProcessGrid(MPI_Comm world_comm, int process_slices,
                                bool be_verbose = false);
//! Construct the global process grid from comm world
//! \param[in] process_slices number of grid slices.
//! \param[in] be_verbose verbosity flag.
void ConstructGlobalProcessGrid(int process_slices, bool be_verbose = false);
//! Get the slice of the current process.
int GetGlobalMySlice();
//! Get the column of the current process.
int GetGlobalMyColumn();
//! Get the row of the current process.
int GetGlobalMyRow();
//! Get the number of process slices.
int GetGlobalNumSlices();
//! Get the number of process columns.
int GetGlobalNumColumns();
//! Get the number of process rows.
int GetGlobalNumRows();
//! Standard destructor
void DestructGlobalProcessGrid();
} // namespace NTPoly
#endif
