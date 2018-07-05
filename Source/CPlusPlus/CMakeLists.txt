################################################################################
set(Csrc
    ChebyshevSolvers.cc
    DensityMatrixSolvers.cc
    EigenBounds.cc
    EigenSolvers.cc
    ExponentialSolvers.cc
    FixedSolversParameters.cc
    GeometryOptimization.cc
    HermiteSolvers.cc
    InverseSolvers.cc
    IterativeSolversParameters.cc
    LinearSolvers.cc
    LoadBalancer.cc
    MatrixMemoryPool.cc
    MinimizerSolvers.cc
    Permutation.cc
    Polynomial.cc
    PMatrixMemoryPool.cc
    ProcessGrid.cc
    PSMatrix.cc
    real_to_complex.c
    RootSolvers.cc
    SignSolvers.cc
    SMatrix.cc
    SolverBase.cc
    SquareRootSolvers.cc
    TrigonometrySolvers.cc
    TripletList.cc
)

set(Chead
    ChebyshevSolvers.h
    DensityMatrixSolvers.h
    EigenBounds.h
    EigenSolvers.h
    ExponentialSolvers.h
    FixedSolversParameters.h
    GeometryOptimization.h
    HermiteSolvers.h
    InverseSolvers.h
    IterativeSolversParameters.h
    LinearSolvers.h
    LoadBalancer.h
    MatrixMemoryPool.h
    MinimizerSolvers.h
    Permutation.h
    PMatrixMemoryPool.h
    Polynomial.h
    ProcessGrid.h
    PSMatrix.h
    RootSolvers.h
    SignSolvers.h
    SMatrix.h
    SolverBase.h
    SquareRootSolvers.h
    TrigonometrySolvers.h
    TripletList.h
    Triplet.h
)

################################################################################
include_directories(${CMAKE_SOURCE_DIR}/Source/C)

add_library(NTPolyCPP STATIC ${Csrc} ${Chead})
target_link_libraries(NTPolyCPP NTPolyWrapper)

set_target_properties(NTPolyCPP PROPERTIES PUBLIC_HEADER "${Chead}")

install(TARGETS NTPolyCPP
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
