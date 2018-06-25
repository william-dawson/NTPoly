################################################################################
set(Csrc
#     ChebyshevSolvers.cc
#     DensityMatrixSolvers.cc
#     EigenBounds.cc
#     EigenSolvers.cc
#     ExponentialSolvers.cc
#     FixedSolversParameters.cc
#     GeometryOptimization.cc
#     HermiteSolvers.cc
#     InverseSolvers.cc
#     IterativeSolversParameters.cc
#     LinearSolvers.cc
#     LoadBalancer.cc
#     MinimizerSolvers.cc
#     Permutation.cc
#     Polynomial.cc
    ProcessGrid.cc
    real_to_complex.c
#     RootSolvers.cc
#     SignSolvers.cc
#     SolverBase.cc
#     SquareRootSolvers.cc
#     TrigonometrySolvers.cc
)

set(Ctypes
  # DistributedSparseMatrix.cc
  # DistributedMatrixMemoryPool.cc
  MatrixMemoryPool.cc
  SparseMatrix.cc
  TripletList.cc
)

set(Chead
    # ChebyshevSolvers.h
    # DensityMatrixSolvers.h
    # EigenBounds.h
    # EigenSolvers.h
    # ExponentialSolvers.h
    # FixedSolversParameters.h
    # GeometryOptimization.h
    # HermiteSolvers.h
    # InverseSolvers.h
    # IterativeSolversParameters.h
    # LinearSolvers.h
    # LoadBalancer.h
    # MinimizerSolvers.h
    # Permutation.h
    # Polynomial.h
    ProcessGrid.h
    # RootSolvers.h
    # SignSolvers.h
    # SolverBase.h
    # SquareRootSolvers.h
    # TrigonometrySolvers.h
)

set(CtypesH
  # DistributedSparseMatrix.h
  # DistributedMatrixMemoryPool.h
  MatrixMemoryPool.h
  SparseMatrix.h
  TripletList.h
  Triplet.h
)

################################################################################
include_directories(${CMAKE_SOURCE_DIR}/Source/C)

add_library(NTPolyCPP STATIC ${Csrc} ${Ctypes} ${Chead} ${CheadH})
target_link_libraries(NTPolyCPP NTPolyWrapper)

set_target_properties(NTPolyCPP PROPERTIES PUBLIC_HEADER "${Chead} ${CheadH}")

install(TARGETS NTPolyCPP
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
