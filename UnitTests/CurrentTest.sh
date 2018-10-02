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
export GEOMH1=@CMAKE_SOURCE_DIR@/UnitTests/Data/F1.mtx
export GEOMO1=@CMAKE_SOURCE_DIR@/UnitTests/Data/S1.mtx
export GEOMO2=@CMAKE_SOURCE_DIR@/UnitTests/Data/S2.mtx
export GEOMD2=@CMAKE_SOURCE_DIR@/UnitTests/Data/D2.mtx
export REALIO=@CMAKE_SOURCE_DIR@/UnitTests/Data/realio.mtx
export CholTest=@CMAKE_SOURCE_DIR@/UnitTests/Data/CholTest.mtx
export HCOMPLEX=@CMAKE_SOURCE_DIR@/UnitTests/Data/complexH.mtx
export SCOMPLEX=@CMAKE_SOURCE_DIR@/UnitTests/Data/complexS.mtx
export DCOMPLEX=@CMAKE_SOURCE_DIR@/UnitTests/Data/complexD.mtx
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

TestName="$(git rev-parse --abbrev-ref HEAD)"
BranchFile="@CMAKE_SOURCE_DIR@/UnitTests/$TestName.sh"

if [ -f $BranchFile ]
then
  source "@CMAKE_SOURCE_DIR@/UnitTests/$TestName.sh"
  @MPIEXEC@ @MPIEXEC_NUMPROC_FLAG@ $PROCESSES @oversubscribe@ \
  @PYTHON_EXECUTABLE@ -m unittest -v $BRANCHTEST
else
  echo "No local testfile ${BranchFile}"
fi
