################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project("Basic Example")
set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

################################################################################
## Packages
find_package(MPI REQUIRED)
find_package(NTPoly REQUIRED)
find_package(NTPolyWrapper REQUIRED)
find_package(NTPolyCPP REQUIRED)
find_package(OpenMP)
find_package(BLAS)
find_package(LAPACK)

################################################################################
## Output Locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

################################################################################
## Fortran Version
add_executable(example_fortran ${CMAKE_SOURCE_DIR}/../PremadeMatrix/main.f90)
target_link_libraries(example_fortran NTPoly::NTPoly ${MPI_Fortran_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${LAPACK_LIBRARIES} 
                      ${BLAS_LIBRARIES})
target_include_directories(example_fortran PRIVATE ${MPI_INCLUDE_PATH})
install(TARGETS example_fortran DESTINATION bin)

## C++ Version
add_executable(example_cpp ${CMAKE_SOURCE_DIR}/../PremadeMatrix/main.cc)
target_link_libraries(example_cpp NTPoly::NTPolyWrapper NTPoly::NTPolyCPP
                      ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES})
install(TARGETS example_cpp DESTINATION bin)
