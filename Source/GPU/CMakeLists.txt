# Library-owned optional backend; no mini-app is required to supply its objects.
set(NTPOLY_CUDA_ARCH sm_100 CACHE STRING "CUDA architecture for the optional cuSOLVERMp backend")
find_program(NTPOLY_EIGEN_NVCC nvcc HINTS $ENV{CUDA_HOME}/bin $ENV{CUDA_PATH}/bin
  /shared/software/hpc_sdk/Linux_aarch64/26.3/cuda/bin)
find_path(NTPOLY_CUSOLVERMP_INCLUDE cusolverMp.h HINTS $ENV{CUSOLVERMP_ROOT}/include
  /shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/include)
find_library(NTPOLY_CUSOLVERMP_LIBRARY cusolverMp HINTS $ENV{CUSOLVERMP_ROOT}/lib64 $ENV{CUSOLVERMP_ROOT}/lib
  /shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64)
find_path(NTPOLY_EIGEN_NCCL_INCLUDE nccl.h HINTS $ENV{NCCL_HOME}/include
  /shared/software/hpc_sdk/Linux_aarch64/26.3/comm_libs/13.1/nccl/include)
find_library(NTPOLY_EIGEN_NCCL_LIBRARY nccl HINTS $ENV{NCCL_HOME}/lib $ENV{NCCL_HOME}/lib64
  /shared/software/hpc_sdk/Linux_aarch64/26.3/comm_libs/13.1/nccl/lib)
if(NOT NTPOLY_EIGEN_NVCC OR NOT NTPOLY_CUSOLVERMP_INCLUDE OR NOT NTPOLY_CUSOLVERMP_LIBRARY
   OR NOT NTPOLY_EIGEN_NCCL_INCLUDE OR NOT NTPOLY_EIGEN_NCCL_LIBRARY)
  message(FATAL_ERROR "NTPOLY_ENABLE_CUSOLVERMP requires nvcc, cuSOLVERMp and NCCL; set CUDA_HOME, CUSOLVERMP_ROOT and NCCL_HOME")
endif()
get_filename_component(NTPOLY_EIGEN_CUDA_BIN ${NTPOLY_EIGEN_NVCC} DIRECTORY)
get_filename_component(NTPOLY_EIGEN_CUDA_ROOT ${NTPOLY_EIGEN_CUDA_BIN}/.. REALPATH)
find_library(NTPOLY_EIGEN_CUDART cudart HINTS ${NTPOLY_EIGEN_CUDA_ROOT}/lib64 REQUIRED)
# NVHPC toolchain files use MPI compiler wrappers; find MPI includes for nvcc.
find_package(MPI REQUIRED COMPONENTS CXX)
set(EIGEN_MPI_FLAGS)
foreach(dir IN LISTS MPI_CXX_INCLUDE_DIRS)
  list(APPEND EIGEN_MPI_FLAGS -I${dir})
endforeach()
set(EIGEN_OBJECT ${CMAKE_CURRENT_BINARY_DIR}/cusolvermp_eigh.o)
add_custom_command(OUTPUT ${EIGEN_OBJECT}
  COMMAND ${NTPOLY_EIGEN_NVCC} -O3 -std=c++17 -ccbin ${CMAKE_CXX_COMPILER} -arch=${NTPOLY_CUDA_ARCH}
    ${EIGEN_MPI_FLAGS} -I${NTPOLY_CUSOLVERMP_INCLUDE} -I${NTPOLY_EIGEN_NCCL_INCLUDE}
    -c ${CMAKE_CURRENT_SOURCE_DIR}/cusolvermp_eigh.cu -o ${EIGEN_OBJECT}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cusolvermp_eigh.cu VERBATIM)
set_source_files_properties(${EIGEN_OBJECT} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE)
add_custom_target(ntpoly_eigen_gpu_object DEPENDS ${EIGEN_OBJECT})
add_dependencies(NTPoly ntpoly_eigen_gpu_object)
target_sources(NTPoly PRIVATE ${EIGEN_OBJECT})
target_link_libraries(NTPoly PUBLIC ${NTPOLY_CUSOLVERMP_LIBRARY} ${NTPOLY_EIGEN_NCCL_LIBRARY}
  ${NTPOLY_EIGEN_CUDART} ${MPI_CXX_LIBRARIES} stdc++)
