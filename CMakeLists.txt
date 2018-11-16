################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project(NTPoly VERSION 2.0)
enable_language(Fortran)

################################################################################
## Packages
if (NOT CMAKE_TOOLCHAIN_FILE)
  message(WARNING "Building without a toolchain file. "
          "If this does not work, please see the example toolchain files in "
          "the Targets directory and set the parameters that match your "
          "system.")
  find_package(MPI REQUIRED)
  find_package(BLAS)
endif()
find_package(SWIG 3.0)

################################################################################
## Testing
if (NOT FORTRAN_ONLY)
  enable_testing()
else()
  message(WARNING "Fortran only! No local testing will be generated.")
endif()

################################################################################
## Output Locations
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scratch)

################################################################################
## Compiler Flags
if (CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_Fortran_FLAGS_DEBUG ${F_TOOLCHAINFLAGS_DEBUG})
  set(CMAKE_CXX_FLAGS_DEBUG ${CXX_TOOLCHAINFLAGS_DEBUG})
  set(CMAKE_Fortran_FLAGS_RELEASE ${F_TOOLCHAINFLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS_RELEASE ${CXX_TOOLCHAINFLAGS_RELEASE})
else()
  find_package(OpenMP)
  set(CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

################################################################################
## Check MPI Features
option(USE_MPIH OFF)
if (USE_MPIH OR NOT ${MPI_Fortran_HAVE_F90_MODULE})
  add_definitions(-DUSE_MPIH)
endif()

option(NOIALLGATHER OFF)
if (NOIALLGATHER)
  add_definitions(-DNOIALLGATHER)
  message(STATUS "IAllgather replacement activated.")
  message(STATUS "Note that this may reduce parallel performance.")
endif()

################################################################################
## Get at the source files
subdirs(Documentation)
subdirs(Source)

if (NOT FORTRAN_ONLY)
  subdirs(UnitTests)
endif()
