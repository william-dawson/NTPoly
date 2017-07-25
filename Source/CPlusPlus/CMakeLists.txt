set(Csrc
    ChebyshevSolvers.cc ChebyshevSolvers.h
    DensityMatrixSolvers.cc DensityMatrixSolvers.h
    DistributedSparseMatrix.cc DistributedSparseMatrix.h
    DistributedMatrixMemoryPool.cc DistributedMatrixMemoryPool.h
    ExponentialSolvers.cc ExponentialSolvers.h
    FixedSolversParameters.cc FixedSolversParameters.h
    InverseSolvers.cc InverseSolvers.h
    IterativeSolversParameters.cc IterativeSolversParameters.h
    LinearSolvers.cc LinearSolvers.h
    LoadBalancer.cc LoadBalancer.h
    MatrixMemoryPool.cc MatrixMemoryPool.h
    MinimizerSolvers.cc MinimizerSolvers.h
    Permutation.cc Permutation.h
    Polynomial.cc Polynomial.h
    ProcessGrid.cc ProcessGrid.h
    RootSolvers.cc RootSolvers.h
    SignSolvers.cc SignSolvers.h
    SolverBase.cc SolverBase.h
    SparseMatrix.cc SparseMatrix.h
    SquareRootSolvers.cc SquareRootSolvers.h
    TrigonometrySolvers.cc TrigonometrySolvers.h
    TripletList.cc TripletList.h
    Triplet.h
)

add_library(NTPolyCPP STATIC ${Csrc})
target_link_libraries(NTPolyCPP NTPolyWrapper)
include_directories(${CMAKE_SOURCE_DIR}/Source/C)
