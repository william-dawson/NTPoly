#!/usr/bin/bash

set -e

################################################################################
## CMake setup
export PYTHONPATH=$PYTHONPATH:@CMAKE_BINARY_DIR@/python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:@CMAKE_BINARY_DIR@/lib
export SCRATCHDIR=@CMAKE_BINARY_DIR@/scratch
export GEOMH1=@CMAKE_SOURCE_DIR@/UnitTests/Data/F1.mtx
export GEOMO1=@CMAKE_SOURCE_DIR@/UnitTests/Data/S1.mtx
export GEOMO2=@CMAKE_SOURCE_DIR@/UnitTests/Data/S2.mtx
export GEOMD2=@CMAKE_SOURCE_DIR@/UnitTests/Data/D2.mtx
export REALIO=@CMAKE_SOURCE_DIR@/UnitTests/Data/realio.mtx
export CholTest=@CMAKE_SOURCE_DIR@/UnitTests/Data/CholTest.mtx
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
  @Python_EXECUTABLE@ -m unittest -v test_matrix
fi

## Matrix Tests
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @oversubscribe@ \
@Python_EXECUTABLE@ -m unittest -v test_solvers
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @oversubscribe@ \
@Python_EXECUTABLE@ -m unittest -v test_chemistry
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @oversubscribe@ \
@Python_EXECUTABLE@ -m unittest -v test_psmatrix
@MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @oversubscribe@ \
@Python_EXECUTABLE@ -m unittest -v test_psmatrixalgebra
