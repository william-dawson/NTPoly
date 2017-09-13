################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project(NTPoly)
set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)
find_package(MPI REQUIRED)

################################################################################
## Testing
enable_testing()

################################################################################
## Output Locations
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scratch)

################################################################################
## Documentation
find_package(Doxygen REQUIRED)
configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in
               ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)

################################################################################
## Toolchain File
set(CMAKE_Fortran_FLAGS ${F_TOOLCHAINFLAGS})
set(CMAKE_CXX_FLAGS ${CXX_TOOLCHAINFLAGS})

################################################################################
## Get at the source files
subdirs(Source)
subdirs(UnitTests)
