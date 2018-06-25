subdirs(Communication)
subdirs(DistributedTypes)
subdirs(Helpers)
subdirs(LocalTypes)

include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Communication)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/DistributedTypes)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/Helpers)
include_directories(${CMAKE_BINARY_DIR}/Source/Fortran/LocalTypes)

# set(Fsrc
#     ChebyshevSolversModule.F90
#     DensityMatrixSolversModule.F90
#     EigenBoundsModule.F90
#     EigenSolversModule.F90
#     ExponentialSolversModule.F90
#     FixedSolversModule.F90
#     GeometryOptimizationModule.F90
#     HermiteSolversModule.F90
#     InverseSolversModule.F90
#     IterativeSolversModule.F90
#     LinearSolversModule.F90
#     LoadBalancerModule.F90
#     MinimizerSolversModule.F90
#     ParameterConverterModule.F90
#     PermutationModule.F90
#     PolynomialSolversModule.F90
#     RootSolversModule.F90
#     SignSolversModule.F90
#     SquareRootSolversModule.F90
#     TrigonometrySolversModule.F90
# )

set(FCommunication
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/MatrixReduceModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Communication/ProcessGridModule.F90
)

# set(FDistributed
#     ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/MatrixMemoryPoolPModule.F90
#     ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/MatrixPModule.F90
#     ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/MatrixPSAlgebraModule.F90
#     ${CMAKE_CURRENT_SOURCE_DIR}/DistributedTypes/MatrixPSModule.F90
# )

set(FHelpers
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/ErrorModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/GemmTasksModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/LoggingModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/MatrixMarketModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Helpers/TimerModule.F90
)

set(FTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/DataTypesModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixDAlgebraModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixDModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixMemoryPoolModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixSAlgebraModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/MatrixSModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/VectorSModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletListModule.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/LocalTypes/TripletModule.F90
)

IF (EIGEN_EXA)
  LIST(APPEND Fsrc EigenExaModule.F90)
  include_directories(${EIGEN_EXA})
ENDIF()

add_library(NTPoly STATIC ${Fsrc} ${FCommunication} ${FDistributed} ${FHelpers}
            ${FTypes})
install(TARGETS NTPoly
        ARCHIVE DESTINATION lib)
