################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.9)
project(NTPoly VERSION 2.7.0 DESCRIPTION 
	"A parallel library for computing the functions of sparse matrices")
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
## Options
option(FORTRAN_ONLY "Build only the Fortran bindings." OFF)
option(NOSWIG "Don't build the swig bindings." OFF)
option(USE_MPIH "Use mpi.h for systems that don't have the mpi module" OFF)
option(NOIALLGATHER
  "A (slower) alternative communication strategy if iallgather not implemented"
  OFF)

## Process these options
if (NOT SWIG_FOUND OR FORTRAN_ONLY)
  set(NOSWIG TRUE)
endif()

################################################################################
## Testing
if (NOT NOSWIG)
  enable_testing()
else()
  message(WARNING
          "Swig disabled. No python bindings or testing will be generated.")
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
if (USE_MPIH OR NOT ${MPI_Fortran_HAVE_F90_MODULE})
  add_definitions(-DUSE_MPIH)
endif()

if (NOIALLGATHER)
  add_definitions(-DNOIALLGATHER)
  message(STATUS "IAllgather replacement activated.")
  message(STATUS "Note that this may reduce parallel performance.")
endif()

################################################################################
## Get at the source files
subdirs(Documentation)
subdirs(Source)

if (NOT NOSWIG)
  subdirs(UnitTests)
endif()

################################################################################
## CMake Packaging
set(CPACK_PACKAGE_VENDOR "RIKEN R-CCS")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES /.git /Build)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.md")
include(CPACK)

