################################################################################
## Basic SETup
cmake_minimum_required (VERSION 3.2)
project(NTPoly VERSION 2.0)
enable_language(Fortran)

################################################################################
## Packages
find_package(MPI REQUIRED)
MESSAGE(STATUS "MPI PATH")
MESSAGE(STATUS ${MPI_INCLUDE_PATH})
find_package(BLAS)
find_package(SWIG 3.0)

################################################################################
## Testing
IF (NOT FORTRAN_ONLY)
  enable_testing()
ELSE()
  MESSAGE(WARNING "Fortran only! No local testing will be generated.")
ENDIF()

################################################################################
## Output Locations
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scratch)

################################################################################
## Compiler Flags From A Toolchain File
option(CMAKE_TOOLCHAIN_FILE OFF)
IF (CMAKE_TOOLCHAIN_FILE)
  SET(CMAKE_Fortran_FLAGS_DEBUG ${F_TOOLCHAINFLAGS_DEBUG})
  SET(CMAKE_CXX_FLAGS_DEBUG ${CXX_TOOLCHAINFLAGS_DEBUG})
  SET(CMAKE_Fortran_FLAGS_RELEASE ${F_TOOLCHAINFLAGS_RELEASE})
  SET(CMAKE_CXX_FLAGS_RELEASE ${CXX_TOOLCHAINFLAGS_RELEASE})
ELSE()
  find_package(OpenMP)
  SET(CMAKE_Fortran_FLAGS_RELEASE
      "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
ENDIF()

## Fill in the default flags if no release type.
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
ENDIF()

################################################################################
## Check MPI Features
option(USE_MPIH OFF)
IF (USE_MPIH OR NOT ${MPI_Fortran_HAVE_F90_MODULE})
  add_definitions(-DUSE_MPIH)
ENDIF()

option(NOIALLGATHER OFF)
IF (NOIALLGATHER)
  add_definitions(-DNOIALLGATHER)
  MESSAGE(STATUS "IAllgather replacement activated.")
  MESSAGE(STATUS "Note that this may reduce parallel performance.")
ENDIF()

################################################################################
## Get at the source files
subdirs(Documentation)
subdirs(Source)

IF (NOT FORTRAN_ONLY)
  subdirs(UnitTests)
ENDIF()
