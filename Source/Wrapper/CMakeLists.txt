subdirs(Communication)
subdirs(DistributedTypes)
subdirs(LocalTypes)

include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/Communication)
include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/DistributedTypes)
include_directories(${CMAKE_BINARY_DIR}/Source/Wrapper/LocalTypes)

set(Wsrc
    ChebyshevSolversModule_wrp.f90
    DensityMatrixSolversModule_wrp.f90
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
    MinimizerSolversModule_wrp.f90
    PermutationModule_wrp.f90
    PolynomialSolversModule_wrp.f90
    RootSolversModule_wrp.f90
    SignSolversModule_wrp.f90
    SquareRootSolversModule_wrp.f90
    TrigonometrySolversModule_wrp.f90
    WrapperModule.f90
)

set(WCommunication
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/ProcessGridModule_wrp.f90
)

set(WDistributed
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PMatrixMemoryPoolModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PSMatrixAlgebraModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/PSMatrixModule_wrp.f90
)

set(WTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixMemoryPoolModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SMatrixAlgebraModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/SMatrixModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletListModule_wrp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletModule_wrp.f90
)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran)

IF (EIGEN_EXA)
  LIST(APPEND Wsrc EigenExaModule_wrp.f90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPolyWrapper STATIC ${Wsrc} ${WCommunication} ${WDistributed}
            ${WTypes})
target_link_libraries(NTPolyWrapper NTPoly)
install(TARGETS NTPolyWrapper ARCHIVE DESTINATION lib)
