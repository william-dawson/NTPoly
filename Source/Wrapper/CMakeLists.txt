set(Wsrc
    ChebyshevSolversModule_wrp.f90
    DensityMatrixSolversModule_wrp.f90
    DistributedSparseMatrixAlgebraModule_wrp.f90
    DistributedSparseMatrixModule_wrp.f90
    DistributedMatrixMemoryPoolModule_wrp.f90
    EigenBoundsModule_wrp.f90
    EigenSolversModule_wrp.f90
    ExponentialSolversModule_wrp.f90
    FixedSolversModule_wrp.f90
    GeometryOptimizationModule_wrp.f90
    HermiteSolversModule_wrp.f90
    IterativeSolversModule_wrp.f90
    InverseSolversModule_wrp.f90
    LinearSolversModule_wrp.f90
    LoadBalancerModule_wrp.f90
    MatrixMemoryPoolModule_wrp.f90
    MinimizerSolversModule_wrp.f90
    PermutationModule_wrp.f90
    PolynomialSolversModule_wrp.f90
    ProcessGridModule_wrp.f90
    RootSolversModule_wrp.f90
    SignSolversModule_wrp.f90
    SparseMatrixAlgebraModule_wrp.f90
    SparseMatrixModule_wrp.f90
    SquareRootSolversModule_wrp.f90
    TrigonometrySolversModule_wrp.f90
    TripletListModule_wrp.f90
    TripletModule_wrp.f90
    WrapperModule.f90
)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran)

IF (EIGEN_EXA)
  LIST(APPEND Wsrc EigenExaModule_wrp.f90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPolyWrapper STATIC ${Wsrc})
target_link_libraries(NTPolyWrapper NTPoly)
install(TARGETS NTPolyWrapper ARCHIVE DESTINATION lib)
