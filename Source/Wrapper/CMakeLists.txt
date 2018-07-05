subdirs(Communication)
subdirs(DistributedTypes)
subdirs(LocalTypes)

include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/Communication)
include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/DistributedTypes)
include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/LocalTypes)

set(Wsrc
    ChebyshevSolversModule_wrp.F90
    DensityMatrixSolversModule_wrp.F90
    EigenBoundsModule_wrp.F90
    EigenSolversModule_wrp.F90
    ExponentialSolversModule_wrp.F90
    GeometryOptimizationModule_wrp.F90
    HermiteSolversModule_wrp.F90
    InverseSolversModule_wrp.F90
    LinearSolversModule_wrp.F90
    LoadBalancerModule_wrp.F90
    MinimizerSolversModule_wrp.F90
    PermutationModule_wrp.F90
    PolynomialSolversModule_wrp.F90
    RootSolversModule_wrp.F90
    SignSolversModule_wrp.F90
    SolverParametersModule_wrp.F90
    SquareRootSolversModule_wrp.F90
    TrigonometrySolversModule_wrp.F90
    WrapperModule.F90
)

set(WCommunication
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/ProcessGridModule_wrp.F90
)

set(WDistributed
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PMatrixMemoryPoolModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PSMatrixAlgebraModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PSMatrixModule_wrp.F90
)

set(WTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixMemoryPoolModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SMatrixAlgebraModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SMatrixModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletListModule_wrp.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletModule_wrp.F90
)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran)

IF (EIGEN_EXA)
  LIST(APPEND Wsrc EigenExaModule_wrp.F90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPolyWrapper STATIC ${Wsrc} ${WCommunication} ${WDistributed}
            ${WTypes})
target_link_libraries(NTPolyWrapper NTPoly)
install(TARGETS NTPolyWrapper ARCHIVE DESTINATION lib)
