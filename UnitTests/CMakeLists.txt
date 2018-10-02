################################################################################
set(TestFiles
    __init__.py
    CurrentTest.sh
    helpers.py
    RunTest.sh
    test_chemistry.py
    test_matrix.py
    test_psmatrix.py
    test_psmatrixalgebra.py
    test_solvers.py
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
configure_file(helpers.py helpers.py COPYONLY)
configure_file(RunTest.sh
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/RunTest.sh @ONLY)
configure_file(CurrentTest.sh
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CurrentTest.sh @ONLY)
configure_file(test_chemistry.py test_chemistry.py COPYONLY)
configure_file(test_matrix.py test_matrix.py COPYONLY)
configure_file(test_psmatrix.py test_psmatrix.py COPYONLY)
configure_file(test_psmatrixalgebra.py test_psmatrixalgebra.py COPYONLY)
configure_file(test_solvers.py test_solvers.py COPYONLY)

################################################################################
find_program(BASH_PROGRAM bash)
set(TestScript ${CMAKE_BINARY_DIR}/bin/RunTest.sh)

add_test(Regression111 ${BASH_PROGRAM} ${TestScript} 1 1 1 1)
add_test(Regression211 ${BASH_PROGRAM} ${TestScript} 2 1 1 2)
add_test(Regression121 ${BASH_PROGRAM} ${TestScript} 1 2 1 2)
add_test(Regression112 ${BASH_PROGRAM} ${TestScript} 1 1 2 2)
add_test(Regression221 ${BASH_PROGRAM} ${TestScript} 2 2 1 4)
add_test(Regression212 ${BASH_PROGRAM} ${TestScript} 2 1 2 4)
add_test(Regression122 ${BASH_PROGRAM} ${TestScript} 1 2 2 4)
add_test(Regression222 ${BASH_PROGRAM} ${TestScript} 2 2 2 8)
add_test(Regression321 ${BASH_PROGRAM} ${TestScript} 3 2 1 6)
add_test(Regression213 ${BASH_PROGRAM} ${TestScript} 2 1 3 6)
add_test(Regression611 ${BASH_PROGRAM} ${TestScript} 6 1 1 6)

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
add_test(CurrentTest311 ${BASH_PROGRAM} ${CurrentTestScript} 3 1 1 3)
add_test(CurrentTest131 ${BASH_PROGRAM} ${CurrentTestScript} 1 3 1 3)
add_test(CurrentTest113 ${BASH_PROGRAM} ${CurrentTestScript} 1 1 3 3)
