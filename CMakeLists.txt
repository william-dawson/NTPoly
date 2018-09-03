################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project(NTPoly VERSION 1.3)
set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

################################################################################
## Packages
find_package(MPI REQUIRED)
find_package(SWIG 3.0)

################################################################################
## Testing
if (NOT FORTRAN_ONLY)
  enable_testing()
else()
  MESSAGE("Fortran only! No local testing will be generated.")
endif()

################################################################################
## Output Locations
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scratch)

################################################################################
## Compiler Flags
set(CMAKE_Fortran_FLAGS_DEBUG ${F_TOOLCHAINFLAGS_DEBUG})
set(CMAKE_CXX_FLAGS_DEBUG ${CXX_TOOLCHAINFLAGS_DEBUG})
set(CMAKE_Fortran_FLAGS_RELEASE ${F_TOOLCHAINFLAGS_RELEASE})
set(CMAKE_CXX_FLAGS_RELEASE ${CXX_TOOLCHAINFLAGS_RELEASE})

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_Fortran_FLAGS ${F_TOOLCHAINFLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS ${CXX_TOOLCHAINFLAGS_RELEASE})
endif()

################################################################################
## Check for EigenExa
IF (EIGEN_EXA)
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -DEIGENEXA")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DEIGENEXA")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DEIGENEXA")

  set(TOOLCHAIN_LIBS "${TOOLCHAIN_LIBS} -L${EIGEN_EXA}")
  set(TOOLCHAIN_LIBS "${TOOLCHAIN_LIBS} -lEigenExa")
ENDIF()

################################################################################
## Get at the source files
subdirs(Documentation)
subdirs(Source)

if (NOT FORTRAN_ONLY)
  subdirs(UnitTests)
endif()
