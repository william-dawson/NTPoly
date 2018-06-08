subdirs(Algebra)
subdirs(Communication)
subdirs(Helpers)
subdirs(Types)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Algebra)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Communication)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Helpers)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Types)

set(Fsrc
    ChebyshevSolversModule.f90
    DensityMatrixSolversModule.f90
    EigenSolversModule.f90
    ExponentialSolversModule.f90
    FixedSolversModule.f90
    GeometryOptimizationModule.f90
    HermiteSolversModule.f90
    InverseSolversModule.f90
    IterativeSolversModule.f90
    LinearSolversModule.f90
    MinimizerSolversModule.f90
    PermutationModule.f90
    PolynomialSolversModule.f90
    RootSolversModule.f90
    SignSolversModule.f90
    SquareRootSolversModule.f90
    TrigonometrySolversModule.f90
)

set(FAlgebra
    ${CMAKE_CURRENT_SOURCE_DIR}/Algebra/DistributedSparseMatrixAlgebraModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Algebra/EigenBoundsModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Algebra/SparseMatrixAlgebraModule.f90
)

set(FCommunication
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixBroadcastModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixGatherModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixMemoryPoolModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixReduceModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixSendRecvModule.f90
)

set(FHelpers
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/GemmTasksModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/LoadBalancerModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/LoggingModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/MatrixMarketModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/TimerModule.f90
)

set(FTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/DataTypesModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/DenseMatrixModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/DistributedMatrixMemoryPoolModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/DistributedSparseMatrixModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/ErrorModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/ParameterConverterModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/ProcessGridModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/SparseMatrixModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/SparseVectorModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/TripletListModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/Types/TripletModule.f90
)

IF (EIGEN_EXA)
  LIST(APPEND Fsrc EigenExaModule.f90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPoly STATIC ${Fsrc} ${FAlgebra} ${FCommunication} ${FHelpers}
            ${FTypes})
install(TARGETS NTPoly
        ARCHIVE DESTINATION lib)
