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
  //! \param[in] be_verbose verbosity flag.
  ProcessGrid(MPI_Comm world_comm, int process_rows, int process_columns,
              int process_slices, bool be_verbose = false);
  //! Construct the process grid from comm world
  //! \param[in] process_rows number of grid rows.
  //! \param[in] process_columns number of grid columns.
  //! \param[in] process_slices number of grid slices.
  //! \param[in] be_verbose verbosity flag.
  ProcessGrid(int process_rows, int process_columns, int process_slices,
              bool be_verbose = false);
  //! Copy constructor.
  //!\param grid to copy from.
  ProcessGrid(const ProcessGrid &old_grid);

public:
  //! Get the slice of the current process.
  int GetMySlice();
  //! Get the column of the current process.
  int GetMyColumn();
  //! Get the row of the current process.
  int GetMyRow();

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
//! Get the slice of the current process.
int GetGlobalMySlice();
//! Get the column of the current process.
int GetGlobalMyColumn();
//! Get the row of the current process.
int GetGlobalMyRow();
//! Standard destructor
void DestructGlobalProcessGrid();
} // namespace NTPoly
#endif
