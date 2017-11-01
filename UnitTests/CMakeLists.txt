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
execute_process(COMMAND ${MPIEXEC} --version
                OUTPUT_VARIABLE mpiversion OUTPUT_STRIP_TRAILING_WHITESPACE)
if (${mpiversion} MATCHES "OpenRTE")
  set(oversubscribe "--oversubscribe")
else()
  set(oversubscribe "")
endif()

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
add_test(PythonTests321 ${BASH_PROGRAM} ${TestScript} 3 2 1 6)
add_test(PythonTests231 ${BASH_PROGRAM} ${TestScript} 2 3 1 6)
add_test(PythonTests213 ${BASH_PROGRAM} ${TestScript} 2 1 3 6)
add_test(PythonTests611 ${BASH_PROGRAM} ${TestScript} 6 1 1 6)

set(CurrentTestScript ${CMAKE_BINARY_DIR}/bin/CurrentTest.sh)
add_test(CurrentTest ${BASH_PROGRAM} ${CurrentTestScript} 1 1 1 1)
add_test(CurrentTest111 ${BASH_PROGRAM} ${CurrentTestScript} 1 1 1 1)
add_test(CurrentTest211 ${BASH_PROGRAM} ${CurrentTestScript} 2 1 1 2)
add_test(CurrentTest121 ${BASH_PROGRAM} ${CurrentTestScript} 1 2 1 2)
add_test(CurrentTest112 ${BASH_PROGRAM} ${CurrentTestScript} 1 1 2 2)
add_test(CurrentTest221 ${BASH_PROGRAM} ${CurrentTestScript} 2 2 1 4)
add_test(CurrentTest212 ${BASH_PROGRAM} ${CurrentTestScript} 2 1 2 4)
add_test(CurrentTest122 ${BASH_PROGRAM} ${CurrentTestScript} 1 2 2 4)
add_test(CurrentTest222 ${BASH_PROGRAM} ${CurrentTestScript} 2 2 2 8)
add_test(CurrentTest321 ${BASH_PROGRAM} ${CurrentTestScript} 3 2 1 6)
add_test(CurrentTest231 ${BASH_PROGRAM} ${CurrentTestScript} 2 3 1 6)
add_test(CurrentTest213 ${BASH_PROGRAM} ${CurrentTestScript} 2 1 3 6)
