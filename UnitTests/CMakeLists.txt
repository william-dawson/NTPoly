################################################################################
set(TestFiles
    __init__.py
    CurrentTest.sh
    Helpers.py
    RunTest.sh
    testChemistry.py
    testDistributedSparseMatrix.py
    testDistributedSparseMatrixAlgebra.py
    testSolvers.py
    testSparseMatrix.py
)

find_package(PythonInterp REQUIRED)

################################################################################
configure_file(__init__.py __init__.py COPYONLY)
configure_file(Helpers.py Helpers.py COPYONLY)
configure_file(RunTest.sh
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/RunTest.sh @ONLY)
configure_file(CurrentTest.sh
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CurrentTest.sh @ONLY)
configure_file(testChemistry.py testChemistry.py COPYONLY)
configure_file(testDistributedSparseMatrix.py
  testDistributedSparseMatrix.py COPYONLY)
configure_file(testDistributedSparseMatrixAlgebra.py
  testDistributedSparseMatrixAlgebra.py COPYONLY)
configure_file(testSolvers.py testSolvers.py COPYONLY)
configure_file(testSparseMatrix.py testSparseMatrix.py COPYONLY)

################################################################################
find_program(BASH_PROGRAM bash)
set(TestScript ${CMAKE_BINARY_DIR}/bin/RunTest.sh)
add_test(PythonTests111 ${BASH_PROGRAM} ${TestScript} 1 1 1 1)
add_test(PythonTests211 ${BASH_PROGRAM} ${TestScript} 2 1 1 2)
add_test(PythonTests121 ${BASH_PROGRAM} ${TestScript} 1 2 1 2)
add_test(PythonTests112 ${BASH_PROGRAM} ${TestScript} 1 1 2 2)
add_test(PythonTests221 ${BASH_PROGRAM} ${TestScript} 2 2 1 4)
add_test(PythonTests212 ${BASH_PROGRAM} ${TestScript} 2 1 2 4)
add_test(PythonTests122 ${BASH_PROGRAM} ${TestScript} 1 2 2 4)
add_test(PythonTests222 ${BASH_PROGRAM} ${TestScript} 2 2 2 8)

set(CurrentTestScript ${CMAKE_BINARY_DIR}/bin/CurrentTest.sh)
add_test(CurrentTest ${BASH_PROGRAM} ${CurrentTestScript} 1 1 1 1)
