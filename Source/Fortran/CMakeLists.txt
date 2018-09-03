################################################################################
set(Fsrc
    ChebyshevSolversModule.F90
    CholeskyModule.F90
    DataTypesModule.F90
    DensityMatrixSolversModule.F90
    DMatrixModule.F90
    EigenBoundsModule.F90
    ErrorModule.F90
    EigenSolversModule.F90
    ExponentialSolversModule.F90
    GemmTasksModule.F90
    GeometryOptimizationModule.F90
    HermiteSolversModule.F90
    InverseSolversModule.F90
    LinearSolversModule.F90
    LoadBalancerModule.F90
    LoggingModule.F90
    MatrixMarketModule.F90
    MatrixMemoryPoolModule.F90
    MatrixReduceModule.F90
    NTMPIModule.F90
    PermutationModule.F90
    PMatrixMemoryPoolModule.F90
    PolynomialSolversModule.F90
    ProcessGridModule.F90
    PSMatrixAlgebraModule.F90
    PSMatrixModule.F90
    RootSolversModule.F90
    SignSolversModule.F90
    SMatrixAlgebraModule.F90
    SMatrixModule.F90
    SolverParametersModule.F90
    SquareRootSolversModule.F90
    SVectorModule.F90
    TimerModule.F90
    TrigonometrySolversModule.F90
    TripletListModule.F90
    TripletModule.F90
)

################################################################################
IF (EIGEN_EXA)
  LIST(APPEND Fsrc EigenExaModule.F90)
  include_directories(${EIGEN_EXA})
ENDIF()

################################################################################
add_library(NTPoly STATIC ${Fsrc})
include(GNUInstallDirs)
install(TARGETS NTPoly
        EXPORT ntpoly-export
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

################################################################################
# For Exporting
install(EXPORT ntpoly-export FILE NTPolyConfig.cmake
        NAMESPACE NTPoly:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NTPoly)
