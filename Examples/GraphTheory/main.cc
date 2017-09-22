////////////////////////////////////////////////////////////////////////////////
// An example based on solving matrices based on premade files.
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include <vector>
using std::vector;
#include <iostream>
// NTPoly Headers
#include "DistributedSparseMatrix.h"
#include "InverseSolvers.h"
#include "IterativeSolversParameters.h"
#include "ProcessGrid.h"
#include "Triplet.h"
#include "TripletList.h"

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // Input Parameters
  string output_file;
  int process_rows, process_columns, process_slices;
  double threshold, convergence_threshold;
  double attenuation;
  int number_of_nodes;
  int extra_connections;

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
    if (key == "--output_file") {
      ss >> output_file;
    } else if (key == "--process_rows") {
      ss >> process_rows;
    } else if (key == "--process_columns") {
      ss >> process_columns;
    } else if (key == "--process_slices") {
      ss >> process_slices;
    } else if (key == "--attenuation") {
      ss >> attenuation;
    } else if (key == "--threshold") {
      ss >> threshold;
    } else if (key == "--number_of_nodes") {
      ss >> number_of_nodes;
    } else if (key == "--extra_connections") {
      ss >> extra_connections;
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
  int number_of_local_nodes = number_of_nodes / total_processors;
  int starting_node = number_of_local_nodes * rank;
  // Handles the edge case
  if (rank == total_processors - 1) {
    number_of_local_nodes = number_of_nodes - rank * number_of_local_nodes;
  }
  vector<int> local_nodes;
  for (int i = 0; i < number_of_local_nodes; ++i) {
    local_nodes.push_back(starting_node + i);
  }
  int ending_node = local_nodes[number_of_local_nodes - 1];

  // Fill The Matrix
  NTPoly::TripletList triplet_list;
  NTPoly::Triplet temp_triplet;

  // First add the connection between each node and itself.
  for (int i = 0; i < number_of_local_nodes; ++i) {
    temp_triplet.index_row = local_nodes[i] + 1;
    temp_triplet.index_column = local_nodes[i] + 1;
    temp_triplet.point_value = 1;
    triplet_list.Append(temp_triplet);
  }

  // Now connections between nearest neighbors.
  for (int i = 0; i < number_of_local_nodes; ++i) {
    temp_triplet.index_row = local_nodes[i] + 1;
    temp_triplet.point_value = 0.1;
    if (local_nodes[i] == 0) {
      // Right Value
      temp_triplet.index_column = local_nodes[i] + 1 + 1;
      triplet_list.Append(temp_triplet);
    } else if (local_nodes[i] == number_of_nodes - 1) {
      // Left Value
      temp_triplet.index_column = local_nodes[i] - 1 + 1;
      triplet_list.Append(temp_triplet);
    } else {
      // Left value
      temp_triplet.index_column = local_nodes[i] - 1 + 1;
      triplet_list.Append(temp_triplet);
      // Right value
      temp_triplet.index_column = local_nodes[i] + 1 + 1;
      triplet_list.Append(temp_triplet);
    }
  }

  // Finally the random extra connections.
  vector<int> extra_scratch(number_of_nodes, 0);
  int counter = 0;
  while (counter < extra_connections) {
    int extra_source_node = rand() % number_of_nodes;
    int extra_destination_node = rand() % number_of_nodes;
    if (extra_scratch[extra_source_node] != 1 &&
        extra_scratch[extra_destination_node] != 1 &&
        extra_source_node != extra_destination_node &&
        extra_source_node != extra_destination_node - 1 &&
        extra_source_node != extra_destination_node + 1) {
      counter = counter + 1;
      extra_scratch[extra_source_node] = 1;
      extra_scratch[extra_destination_node] = 1;
      if (extra_source_node >= starting_node &&
          extra_source_node <= ending_node) {
        temp_triplet.index_row = extra_source_node + 1;
        temp_triplet.index_column = extra_destination_node + 1;
        temp_triplet.point_value = 0.1;
        triplet_list.Append(temp_triplet);
      } else if (extra_destination_node >= starting_node &&
                 extra_destination_node <= ending_node) {
        temp_triplet.index_row = extra_destination_node + 1;
        temp_triplet.index_column = extra_source_node + 1;
        temp_triplet.point_value = 0.1;
        triplet_list.Append(temp_triplet);
      }
    }
  }

  // Finally build the matrix
  NTPoly::DistributedSparseMatrix NetworkMat(number_of_nodes);
  NetworkMat.FillFromTripletList(triplet_list);

  // Solve
  NTPoly::DistributedSparseMatrix ResMat(number_of_nodes);
  ResMat.FillIdentity();
  ResMat.Increment(NetworkMat, -1.0 * attenuation);

  NTPoly::DistributedSparseMatrix ResultMat(number_of_nodes);
  NTPoly::InverseSolvers::Invert(ResMat, ResultMat, solver_parameters);

  // Print the density matrix to file.
  ResultMat.WriteToMatrixMarket(output_file);

  // Cleanup
  MPI_Finalize();
  return 0;
}
