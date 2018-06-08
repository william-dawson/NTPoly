subdirs(Communication)
subdirs(DistributedTypes)
subdirs(Helpers)
subdirs(LocalTypes)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Communication)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/DistributedTypes)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Helpers)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/LocalTypes)

set(Fsrc
    ChebyshevSolversModule.f90
    DensityMatrixSolversModule.f90
    EigenBoundsModule.f90
    EigenSolversModule.f90
    ExponentialSolversModule.f90
    FixedSolversModule.f90
    GeometryOptimizationModule.f90
    HermiteSolversModule.f90
    InverseSolversModule.f90
    IterativeSolversModule.f90
    LinearSolversModule.f90
    MinimizerSolversModule.f90
    ParameterConverterModule.f90
    PermutationModule.f90
    PolynomialSolversModule.f90
    RootSolversModule.f90
    SignSolversModule.f90
    SquareRootSolversModule.f90
    TrigonometrySolversModule.f90
)

set(FCommunication
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixBroadcastModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixGatherModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixReduceModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixSendRecvModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/ProcessGridModule.f90
)

set(FDistributed
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/DistributedMatrixMemoryPoolModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/DistributedSparseMatrixAlgebraModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/DistributedSparseMatrixModule.f90
)

set(FHelpers
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/ErrorModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/GemmTasksModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/LoadBalancerModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/LoggingModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/MatrixMarketModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/TimerModule.f90
)

set(FTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/DataTypesModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/DenseMatrixModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixMemoryPoolModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SparseMatrixAlgebraModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SparseMatrixModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SparseVectorModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletListModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletModule.f90
)

IF (EIGEN_EXA)
  LIST(APPEND Fsrc EigenExaModule.f90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPoly STATIC ${Fsrc} ${FCommunication} ${FDistributed} ${FHelpers}
            ${FTypes})
install(TARGETS NTPoly
        ARCHIVE DESTINATION lib)
