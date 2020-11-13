////////////////////////////////////////////////////////////////////////////////
// Study the exponential of a nonsymmetric graph using the technique of
// Guo, Krystal, and Bojan Mohar. "Hermitian adjacency matrix of digraphs and
// mixed graphs." Journal of Graph Theory 85, no. 1 (2017): 217-248.
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include "ExponentialSolvers.h"
#include "Logging.h"
#include "PSMatrix.h"
#include "ProcessGrid.h"
#include "SolverParameters.h"
#include "Triplet.h"
#include "TripletList.h"
#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////
void ConstructGuoMatrix(const NTPoly::Matrix_ps &InMat,
                        NTPoly::Matrix_ps &OutMat) {
  // First Symmetrize The Input Matrix.
  NTPoly::TripletList_r tlist, stlist;
  InMat.GetTripletList(tlist);
  NTPoly::Triplet_r temp, temp2;

  for (int i = 0; i < tlist.GetSize(); ++i) {
    temp = tlist.GetTripletAt(i);
    stlist.Append(temp);
    if (temp.index_row != temp.index_column) {
      temp2.index_row = temp.index_column;
      temp2.index_column = temp.index_row;
      temp2.point_value = temp.point_value;
      stlist.Append(temp2);
    }
  }
  NTPoly::Matrix_ps SMat(InMat.GetActualDimension());
  SMat.FillFromTripletList(stlist);

  // Construct The Guide Matrix.
  NTPoly::Matrix_ps Guide(SMat);
  Guide.Increment(InMat, -1.0);

  // Now iterate over the entries in the guide matrix.
  NTPoly::Triplet_c temp_c;
  NTPoly::TripletList_c clist;
  Guide.GetTripletList(tlist);
  for (int i = 0; i < tlist.GetSize(); ++i) {
    temp = tlist.GetTripletAt(i);
    temp_c.index_row = temp.index_row;
    temp_c.index_column = temp.index_column;
    temp_c.point_value = std::complex<double>(0.0, 1.0);
    clist.Append(temp_c);
  }
  NTPoly::Matrix_ps CMatrix(SMat.GetActualDimension());
  CMatrix.FillFromTripletList(clist);

  // Now add it all together.
  OutMat.Transpose(CMatrix);
  OutMat.Conjugate();
  OutMat.Increment(CMatrix);
  OutMat.Increment(SMat);
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // Input Parameters
  string input_file, exponential_file;
  int process_rows, process_columns, process_slices;
  double threshold;

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
    if (key == "--input_file") {
      ss >> input_file;
    } else if (key == "--exponential_file") {
      ss >> exponential_file;
    } else if (key == "--process_rows") {
      ss >> process_rows;
    } else if (key == "--process_columns") {
      ss >> process_columns;
    } else if (key == "--process_slices") {
      ss >> process_slices;
    } else if (key == "--threshold") {
      ss >> threshold;
    }
  }

  // Setup the process grid.
  NTPoly::ConstructGlobalProcessGrid(MPI_COMM_WORLD, process_rows,
                                     process_columns, process_slices, true);
  if (NTPoly::GetGlobalIsRoot()) {
     NTPoly::ActivateLogger();
  }

  // Compute the Hermitian Matrix
  NTPoly::Matrix_ps InMat(input_file);
  NTPoly::Matrix_ps GMat(InMat.GetActualDimension());
  ConstructGuoMatrix(InMat, GMat);
  GMat.Scale(0.5);

  // Compute The Exponential
  NTPoly::SolverParameters param;
  param.SetThreshold(threshold);
  NTPoly::Matrix_ps OMat(GMat.GetActualDimension());
  NTPoly::ExponentialSolvers::ComputeExponential(GMat, OMat, param);

  // Write To File
  OMat.WriteToMatrixMarket(exponential_file);

  // Cleanup
  if (NTPoly::GetGlobalIsRoot()) {
     NTPoly::DeactivateLogger();
  }
  NTPoly::DestructGlobalProcessGrid();
  MPI_Finalize();
  return 0;
}
