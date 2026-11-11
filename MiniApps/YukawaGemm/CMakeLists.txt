################################################################################
# YukawaGemm mini-app: benchmarks NTPoly's sparse GEMM kernel.
if(NTPOLY_ENABLE_GPU_TRS2)
  add_executable(yukawa_gemm yukawa_gemm.f90)
else()
  add_executable(yukawa_gemm solver_csr.f90 solver_nccl.f90 yukawa_gemm.f90)
endif()
target_include_directories(yukawa_gemm PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(yukawa_gemm PRIVATE NTPoly ${TOOLCHAIN_LIBS})

add_executable(trs2_reference trs2_reference.F90)
target_include_directories(trs2_reference PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(trs2_reference PRIVATE NTPoly ${TOOLCHAIN_LIBS})
add_executable(identity_test tests/test_identity.f90)
target_include_directories(identity_test PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(identity_test PRIVATE NTPoly ${TOOLCHAIN_LIBS})
option(NTPOLY_REFERENCE_OPENACC "Bind reference MPI ranks to OpenACC GPUs" OFF)
if(NTPOLY_REFERENCE_OPENACC)
    target_compile_definitions(trs2_reference PRIVATE REFERENCE_OPENACC)
endif()
option(YUKAWA_CUDA_PROFILE "Enable verbose CUDA events and phase logging" ON)
set(CUSPARSE_PROFILE_FLAG)
if(YUKAWA_CUDA_PROFILE)
    set(CUSPARSE_PROFILE_FLAG -DCUSPARSE_PROFILE)
endif()

# cuSPARSE hybrid backend (CUDA + cuSPARSE library).
# Track CUDA source changes as a normal build dependency.
set(CUSPARSE_ENABLED FALSE)

find_program(NVCC_BIN nvcc
    HINTS $ENV{CUDA_HOME}/bin $ENV{CUDA_PATH}/bin
    PATHS /shared/software/hpc_sdk/Linux_aarch64/26.3/cuda/bin
    NO_DEFAULT_PATH)

if (NVCC_BIN)
    get_filename_component(CUDA_BIN_DIR ${NVCC_BIN} DIRECTORY)
    get_filename_component(CUDA_HOME ${CUDA_BIN_DIR}/.. REALPATH)

    find_library(CUSPARSE_LIB cusparse
        HINTS ${CUDA_HOME}/../math_libs/lib64 ${CUDA_HOME}/lib64
        NO_DEFAULT_PATH)

    if (CUSPARSE_LIB)
        find_path(NCCL_INCLUDE_DIR nccl.h
            HINTS $ENV{NCCL_HOME}/include
            PATHS /shared/software/hpc_sdk/Linux_aarch64/26.3/comm_libs/13.1/nccl/include)
        find_library(NCCL_LIBRARY nccl
            HINTS $ENV{NCCL_HOME}/lib $ENV{NCCL_HOME}/lib64
            PATHS /shared/software/hpc_sdk/Linux_aarch64/26.3/comm_libs/13.1/nccl/lib)
        set(NCCL_COMPILE_FLAGS)
        if(NCCL_INCLUDE_DIR AND NCCL_LIBRARY)
            set(NCCL_COMPILE_FLAGS -DYUKAWA_HAVE_NCCL -I${NCCL_INCLUDE_DIR})
            target_link_libraries(yukawa_gemm PRIVATE ${NCCL_LIBRARY})
            message(STATUS "NCCL solver backend: ENABLED")
        else()
            message(STATUS "NCCL not found; --nccl will report unavailable")
        endif()
        set(NCCL_OBJ ${CMAKE_CURRENT_BINARY_DIR}/nccl_panels.o)
        add_custom_command(OUTPUT ${NCCL_OBJ}
            COMMAND ${NVCC_BIN} -O2 -std=c++17 -ccbin ${CMAKE_CXX_COMPILER}
                -arch=sm_100 ${NCCL_COMPILE_FLAGS}
                -c ${CMAKE_CURRENT_SOURCE_DIR}/nccl_panels.cu -o ${NCCL_OBJ}
                -I${CUDA_HOME}/include
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/nccl_panels.cu
            VERBATIM)
        set_source_files_properties(${NCCL_OBJ} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE)
        if(NTPOLY_ENABLE_GPU_TRS2)
            target_sources(NTPoly PRIVATE ${NCCL_OBJ})
        else()
            target_sources(yukawa_gemm PRIVATE ${NCCL_OBJ})
        endif()
        if(NCCL_INCLUDE_DIR AND NCCL_LIBRARY)
            add_executable(nccl_transport_test tests/nccl_transport_test.cpp ${NCCL_OBJ})
            target_include_directories(nccl_transport_test PRIVATE ${CUDA_HOME}/include)
            target_link_libraries(nccl_transport_test PRIVATE ${NCCL_LIBRARY} ${CUDA_HOME}/lib64/libcudart.so)
        endif()
        set(CUSPARSE_OBJ ${CMAKE_CURRENT_BINARY_DIR}/cusparse_gemm.o)
        add_custom_command(OUTPUT ${CUSPARSE_OBJ}
            COMMAND ${NVCC_BIN} -O2 -std=c++17 -ccbin ${CMAKE_CXX_COMPILER}
                    -arch=sm_100
                    ${CUSPARSE_PROFILE_FLAG}
                    -c ${CMAKE_CURRENT_SOURCE_DIR}/cusparse_gemm.cu
                    -o ${CUSPARSE_OBJ}
                    -I${CUDA_HOME}/include
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cusparse_gemm.cu
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Compiling YukawaGemm CUDA backend"
            VERBATIM)
        set_source_files_properties(${CUSPARSE_OBJ} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE)
        if(NTPOLY_ENABLE_GPU_TRS2)
            target_sources(NTPoly PRIVATE ${CUSPARSE_OBJ})
        else()
            target_sources(yukawa_gemm PRIVATE ${CUSPARSE_OBJ})
        endif()
        set(RESIDENT_OBJ ${CMAKE_CURRENT_BINARY_DIR}/resident_panels.o)
        add_custom_command(OUTPUT ${RESIDENT_OBJ}
            COMMAND ${NVCC_BIN} -O2 -std=c++17 -ccbin ${CMAKE_CXX_COMPILER}
                -arch=sm_100 -c ${CMAKE_CURRENT_SOURCE_DIR}/resident_panels.cu
                -o ${RESIDENT_OBJ} -I${CUDA_HOME}/include
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/resident_panels.cu VERBATIM)
        set_source_files_properties(${RESIDENT_OBJ} PROPERTIES GENERATED TRUE EXTERNAL_OBJECT TRUE)
        if(NTPOLY_ENABLE_GPU_TRS2)
            target_sources(NTPoly PRIVATE ${RESIDENT_OBJ})
        else()
            target_sources(yukawa_gemm PRIVATE ${RESIDENT_OBJ})
        endif()
        if(NCCL_INCLUDE_DIR AND NCCL_LIBRARY)
            add_executable(resident_panels_test tests/resident_panels_test.cpp
                ${RESIDENT_OBJ} ${CUSPARSE_OBJ} ${NCCL_OBJ})
            target_include_directories(resident_panels_test PRIVATE ${CUDA_HOME}/include)
            target_link_libraries(resident_panels_test PRIVATE ${NCCL_LIBRARY}
                ${CUSPARSE_LIB} ${CUDA_HOME}/lib64/libcudart.so)
        endif()
        add_executable(device_partials_test tests/device_partials_test.cpp ${CUSPARSE_OBJ})
        target_link_libraries(device_partials_test PRIVATE ${CUSPARSE_LIB}
            ${CUDA_HOME}/lib64/libcudart.so)
        target_link_libraries(yukawa_gemm PRIVATE ${CUSPARSE_LIB})
        target_link_libraries(yukawa_gemm PRIVATE ${CUDA_HOME}/lib64/libcudart.so)
        target_link_libraries(yukawa_gemm PRIVATE stdc++)
        target_include_directories(yukawa_gemm PRIVATE ${CUDA_HOME}/include)
        set(CUSPARSE_ENABLED TRUE)
        message(STATUS "cuSPARSE hybrid backend: ENABLED")
    else()
        message(STATUS "cuSPARSE library not found — cuSPARSE backend disabled")
    endif()
else()
    message(STATUS "nvcc not found — cuSPARSE backend disabled")
endif()

if (NOT CUSPARSE_ENABLED)
    target_compile_definitions(yukawa_gemm PRIVATE NO_CUSPARSE)
endif()

# Self-contained Yukawa matrix generator (no MPI, no external deps).
add_executable(generate_yukawa generate_yukawa.f90)

if(NTPOLY_ENABLE_GPU_TRS2)
    if(NOT CUSPARSE_ENABLED OR NOT NCCL_INCLUDE_DIR OR NOT NCCL_LIBRARY)
        message(FATAL_ERROR "NTPOLY_ENABLE_GPU_TRS2 requires CUDA, cuSPARSE and NCCL")
    endif()
    add_custom_target(ntpoly_gpu_objects DEPENDS ${CUSPARSE_OBJ} ${NCCL_OBJ} ${RESIDENT_OBJ})
    add_dependencies(NTPoly ntpoly_gpu_objects)
    target_link_libraries(NTPoly PUBLIC ${CUSPARSE_LIB} ${NCCL_LIBRARY}
        ${CUDA_HOME}/lib64/libcudart.so stdc++)
endif()
if(NTPOLY_ENABLE_GPU_TRS2)
    add_executable(native_gpu_test tests/native_gpu_test.F90)
    target_include_directories(native_gpu_test PRIVATE ${CMAKE_BINARY_DIR}/include)
    target_link_libraries(native_gpu_test PRIVATE NTPoly ${TOOLCHAIN_LIBS})
endif()

if(NTPOLY_ENABLE_GPU_TRS2)
    add_executable(gpu_planning_test tests/test_gpu_planning.f90)
    target_include_directories(gpu_planning_test PRIVATE ${CMAKE_BINARY_DIR}/include)
    target_link_libraries(gpu_planning_test PRIVATE NTPoly ${TOOLCHAIN_LIBS})
endif()
