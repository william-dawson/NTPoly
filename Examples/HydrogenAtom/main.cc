////////////////////////////////////////////////////////////////////////////////
// An example based on solving matrices based on a 1D hydrogen molecule.
#include <math.h>
#include <mpi.h>
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include <vector>
using std::vector;
#include <iostream>
// NTPoly Headers
#include "DensityMatrixSolvers.h"
#include "DistributedSparseMatrix.h"
#include "IterativeSolversParameters.h"
#include "ProcessGrid.h"
#include "Triplet.h"
#include "TripletList.h"
// ETC
const double x_start = -6.28;
const double x_end = 6.28;

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // Input Parameters
  string density_file_out;
  int process_rows, process_columns, process_slices;
  double threshold, convergence_threshold;
  int number_of_electrons;
  int grid_points;

  // Setup MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  int total_processors;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processors);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Process The Input Parameters
  string key, value;
  for (int i = 1; i < argc; i += 2) {
    key = string(argv[i]);
    value = string(argv[i + 1]);
    stringstream ss;
    ss << value;
    if (key == "--density") {
      ss >> density_file_out;
    } else if (key == "--process_rows") {
      ss >> process_rows;
    } else if (key == "--process_columns") {
      ss >> process_columns;
    } else if (key == "--process_slices") {
      ss >> process_slices;
    } else if (key == "--grid_points") {
      ss >> grid_points;
    } else if (key == "--number_of_electrons") {
      ss >> number_of_electrons;
    } else if (key == "--threshold") {
      ss >> threshold;
    } else if (key == "--convergence_threshold") {
      ss >> convergence_threshold;
    }
  }

  // Setup the process grid.
  NTPoly::ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns,
                               process_slices, true);
  // Set Up The Solver Parameters.
  NTPoly::IterativeSolverParameters solver_parameters;
  solver_parameters.SetConvergeDiff(convergence_threshold);
  solver_parameters.SetThreshold(threshold);
  solver_parameters.SetVerbosity(true);

  // Divide The Work Amongst Processors.
  int local_grid_points = grid_points / total_processors;
  int start_row = local_grid_points * rank;
  // Handle the edge case
  if (rank == total_processors - 1) {
    local_grid_points = grid_points - rank * local_grid_points;
  }
  vector<double> local_rows(local_grid_points);
  for (int i = 0; i < local_grid_points; ++i) {
    local_rows[i] = start_row + i;
  }

  // Construct A Linear Space.
  vector<double> x_values(local_grid_points);
  double grid_spacing = (x_end - x_start) / (grid_points - 1);
  double local_x_start = x_start + start_row * grid_spacing;
  for (int i = 0; i < local_grid_points; ++i) {
    x_values[i] = local_x_start + i * grid_spacing;
  }

  // Construct The Kinetic Energy Operator.
  NTPoly::TripletList triplet_list;
  NTPoly::Triplet temp_value;
  for (int counter = 0; counter < local_grid_points; ++counter) {
    temp_value.index_row = start_row + counter + 1;
    // Stencil Point 1
    if (temp_value.index_row > 2) {
      temp_value.index_column = temp_value.index_row - 2;
      temp_value.point_value =
          (-0.5) * (-1.0 / (12.0 * grid_spacing * grid_spacing));
      triplet_list.Append(temp_value);
    }
    // Stencil Point 2
    if (temp_value.index_row > 1) {
      temp_value.index_column = temp_value.index_row - 1;
      temp_value.point_value =
          (-0.5) * (16.0 / (12.0 * grid_spacing * grid_spacing));
      triplet_list.Append(temp_value);
    }
    // Stencil Point 3
    temp_value.index_column = temp_value.index_row;
    temp_value.point_value =
        (-0.5) * (-30.0 / (12.0 * grid_spacing * grid_spacing));
    triplet_list.Append(temp_value);
    // Stencil Point 4
    if (temp_value.index_row + 1 < grid_points) {
      temp_value.index_column = temp_value.index_row + 1;
      temp_value.point_value =
          (-0.5) * (16.0 / (12.0 * grid_spacing * grid_spacing));
      triplet_list.Append(temp_value);
    }
    // Stencil Point 5
    if (temp_value.index_row + 2 < grid_points) {
      temp_value.index_column = temp_value.index_row + 2;
      temp_value.point_value =
          (-0.5) * (-1.0 / (12.0 * grid_spacing * grid_spacing));
      triplet_list.Append(temp_value);
    }
  }
  NTPoly::DistributedSparseMatrix KineticEnergy(grid_points);
  KineticEnergy.FillFromTripletList(triplet_list);

  // Construct The Potential Energy Operator.
  NTPoly::TripletList potential_triplet_list;
  for (int i = 0; i < local_grid_points; ++i) {
    temp_value.index_row = start_row + i + 1;
    temp_value.index_column = start_row + i + 1;
    temp_value.point_value = -1.0 / fabs(x_values[i]);
    potential_triplet_list.Append(temp_value);
  }
  NTPoly::DistributedSparseMatrix PotentialEnergy(grid_points);
  PotentialEnergy.FillFromTripletList(potential_triplet_list);

  // Construct The Full Hamiltonian.
  NTPoly::DistributedSparseMatrix Hamiltonian(KineticEnergy);
  Hamiltonian.Increment(PotentialEnergy);

  // Overlap Matrix is just the identity.
  NTPoly::DistributedSparseMatrix Identity(grid_points);
  Identity.FillIdentity();

  // Call the solver routine.
  NTPoly::DistributedSparseMatrix Density(grid_points);
  double chemical_potential;
  NTPoly::DensityMatrixSolvers::TRS2(Hamiltonian, Identity, 2, Density,
                                     chemical_potential, solver_parameters);

  // Print the density matrix to file.
  Density.WriteToMatrixMarket(density_file_out);

  // Cleanup
  NTPoly::DestructProcessGrid();
  MPI_Finalize();
  return 0;
}
