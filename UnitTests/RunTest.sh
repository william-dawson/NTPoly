#!/usr/bin/bash

set -e

################################################################################
## CMake setup
export PYTHONPATH=$PYTHONPATH:@CMAKE_BINARY_DIR@/python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:@CMAKE_BINARY_DIR@/lib
export SCRATCHDIR=@CMAKE_BINARY_DIR@/scratch
export HAMILTONIAN=@CMAKE_SOURCE_DIR@/UnitTests/Data/Hamiltonian.mtx
export OVERLAP=@CMAKE_SOURCE_DIR@/UnitTests/Data/Overlap.mtx
export DENSITY=@CMAKE_SOURCE_DIR@/UnitTests/Data/Density.mtx
cd @CMAKE_BINARY_DIR@/UnitTests

# Get Parameters
if [ "$#" -ne 4 ]
then
  echo "Illegal number of parameters"
  exit
fi
export PROCESS_COLUMNS="$1"
export PROCESS_ROWS="$2"
export PROCESS_SLICES="$3"
export PROCESSES="$4"

## Local Tests
if [ $PROCESSES == "1" ]
then
  @PYTHON_EXECUTABLE@ -m unittest -v testSparseMatrix.TestLocalMatrix
fi

## MPI Tests
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @PYTHON_EXECUTABLE@ \
-m unittest -v testDistributedSparseMatrix.TestDistributedMatrix
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @PYTHON_EXECUTABLE@ \
-m unittest -v testDistributedSparseMatrixAlgebra.TestDistributedMatrixAlgebra
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @PYTHON_EXECUTABLE@ \
-m unittest -v testSolvers.TestSolvers
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @PYTHON_EXECUTABLE@ \
-m unittest -v testChemistry.TestChemistry
