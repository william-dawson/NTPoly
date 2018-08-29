////////////////////////////////////////////////////////////////////////////////
// An example based on solving matrices based on premade files.
#include <mpi.h>
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
// NTPoly Headers
#include "DensityMatrixSolvers.h"
#include "PSMatrix.h"
#include "Permutation.h"
#include "ProcessGrid.h"
#include "SolverParameters.h"
#include "SquareRootSolvers.h"

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // Input Parameters
  string hamiltonian_file;
  string overlap_file;
  string density_file_out;
  int process_rows, process_columns, process_slices;
  double threshold;
  double converge_overlap, converge_density;
  int number_of_electrons;

  // Setup MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

  // Process The Input Parameters
  string key, value;
  for (int i = 1; i < argc; i += 2) {
    key = string(argv[i]);
    value = string(argv[i + 1]);
    stringstream ss;
    ss << value;
    if (key == "--hamiltonian") {
      ss >> hamiltonian_file;
    } else if (key == "--overlap") {
      ss >> overlap_file;
    } else if (key == "--density") {
      ss >> density_file_out;
    } else if (key == "--process_rows") {
      ss >> process_rows;
    } else if (key == "--process_columns") {
      ss >> process_columns;
    } else if (key == "--process_slices") {
      ss >> process_slices;
    } else if (key == "--number_of_electrons") {
      ss >> number_of_electrons;
    } else if (key == "--threshold") {
      ss >> threshold;
    } else if (key == "--converge_overlap") {
      ss >> converge_overlap;
    } else if (key == "--converge_density") {
      ss >> converge_density;
    }
  }

  // Setup the process grid.
  NTPoly::ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns,
                               process_slices, true);

  // Read in the matrices from file.
  NTPoly::Matrix_ps Hamiltonian(hamiltonian_file);
  NTPoly::Matrix_ps Overlap(overlap_file);
  NTPoly::Matrix_ps ISQOverlap(Hamiltonian.GetActualDimension());
  NTPoly::Matrix_ps Density(Hamiltonian.GetActualDimension());

  // Set Up The Solver Parameters.
  NTPoly::Permutation permutation(Hamiltonian.GetLogicalDimension());
  permutation.SetRandomPermutation();
  NTPoly::SolverParameters solver_parameters;
  solver_parameters.SetConvergeDiff(converge_overlap);
  solver_parameters.SetThreshold(threshold);
  solver_parameters.SetLoadBalance(permutation);
  solver_parameters.SetVerbosity(true);

  // Call the solver routines.
  NTPoly::SquareRootSolvers::InverseSquareRoot(Overlap, ISQOverlap,
                                               solver_parameters);

  // Change the solver variable for computing the density matrix.
  solver_parameters.SetConvergeDiff(converge_density);

  // Compute the density matrix.
  double chemical_potential, energy;
  NTPoly::DensityMatrixSolvers::TRS2(Hamiltonian, ISQOverlap,
                                     number_of_electrons, Density, energy,
                                     chemical_potential, solver_parameters);

  // Print the density matrix to file.
  Density.WriteToMatrixMarket(density_file_out);

  // Cleanup
  NTPoly::DestructProcessGrid();
  MPI_Finalize();
  return 0;
}
