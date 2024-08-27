////////////////////////////////////////////////////////////////////////////////
// An example based on solving matrices based on premade files.
#include <mpi.h>
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
// NTPoly Headers
#include "Logging.h"
#include "MatrixMapper.h"
#include "PSMatrix.h"
#include "ProcessGrid.h"

// This is the function we will map on to the matrix.
class TestOperation : public NTPoly::RealOperation {
public:
  // Notice you overload the () operator.
  bool operator()() {
    // This object contains a triplet called data for you to modify.
    if (this->data.index_row >= this->data.index_column) {
      this->data.point_value *= 2;
      return true;
    }
    return false;
  }
};

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // Input Parameters
  string input_file;
  string output_file;
  int process_slices;

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
    if (key == "--input_matrix") {
      ss >> input_file;
    } else if (key == "--output_matrix") {
      ss >> output_file;
    } else if (key == "--process_slices") {
      ss >> process_slices;
    }
  }

  // Setup the process grid.
  NTPoly::ConstructGlobalProcessGrid(MPI_COMM_WORLD, process_slices);
  if (NTPoly::GetGlobalIsRoot()) {
    NTPoly::ActivateLogger();
  }

  // Read in the matrices from file.
  NTPoly::Matrix_ps Input(input_file);
  NTPoly::Matrix_ps Output(Input.GetActualDimension());

  // Map
  TestOperation op;
  NTPoly::MatrixMapper::Map(Input, Output, &op);

  // Print the density matrix to file.
  Output.WriteToMatrixMarket(output_file);

  // Cleanup
  if (NTPoly::GetGlobalIsRoot()) {
    NTPoly::DeactivateLogger();
  }
  NTPoly::DestructGlobalProcessGrid();
  MPI_Finalize();
  return 0;
}
