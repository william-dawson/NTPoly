################################################################################
set(Fsrc
    AnalysisModule.F90
    ChebyshevSolversModule.F90
    CholeskySolversModule.F90
    DataTypesModule.F90
    DensityMatrixSolversModule.F90
    DMatrixModule.F90
    EigenBoundsModule.F90
    EigenExaModule.F90
    EigenSolversModule.F90
    ErrorModule.F90
    ExponentialSolversModule.F90
    FermiOperatorModule.F90
    GemmTasksModule.F90
    GeometryOptimizationModule.F90
    HermiteSolversModule.F90
    InverseSolversModule.F90
    LinearSolversModule.F90
    LoadBalancerModule.F90
    LoggingModule.F90
    MatrixConversionModule.F90
    MatrixMarketModule.F90
    MatrixMapsModule.F90
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
    SingularValueSolversModule.F90
    SolverParametersModule.F90
    SquareRootSolversModule.F90
    SVectorModule.F90
    TimerModule.F90
    TrigonometrySolversModule.F90
    TripletListModule.F90
    TripletModule.F90
)

################################################################################
add_library(NTPoly ${Fsrc})
if (NOT NOSWIG)
  set_target_properties(NTPoly PROPERTIES POSITION_INDEPENDENT_CODE True)
endif()
target_link_libraries(NTPoly ${MPI_Fortran_LIBRARIES} ${EigenSolver_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${LAPACK_LIBRARIES} 
                      ${BLAS_LIBRARIES})
target_include_directories(NTPoly PUBLIC ${EigenSolver_INCLUDE_DIRS}
                           ${MPI_Fortran_INCLUDE_PATH})

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

################################################################################
# pkg-config
configure_file(NTPoly.pc.in NTPoly.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/Source/Fortran/NTPoly.pc 
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

