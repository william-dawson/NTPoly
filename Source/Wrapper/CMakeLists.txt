################################################################################
set(Wsrc
    ChebyshevSolversModule_wrp.F90
    DensityMatrixSolversModule_wrp.F90
    EigenBoundsModule_wrp.F90
    ExponentialSolversModule_wrp.F90
    FermiOperatorExpansionModule_wrp.F90
    GeometryOptimizationModule_wrp.F90
    HermiteSolversModule_wrp.F90
    InverseSolversModule_wrp.F90
    LinearSolversModule_wrp.F90
    LoadBalancerModule_wrp.F90
    MatrixMemoryPoolModule_wrp.F90
    PermutationModule_wrp.F90
    PMatrixMemoryPoolModule_wrp.F90
    PolynomialSolversModule_wrp.F90
    ProcessGridModule_wrp.F90
    PSMatrixAlgebraModule_wrp.F90
    PSMatrixModule_wrp.F90
    RootSolversModule_wrp.F90
    SignSolversModule_wrp.F90
    SMatrixAlgebraModule_wrp.F90
    SMatrixModule_wrp.F90
    SolverParametersModule_wrp.F90
    SquareRootSolversModule_wrp.F90
    TrigonometrySolversModule_wrp.F90
    TripletListModule_wrp.F90
    WrapperModule.F90
)

################################################################################
add_library(NTPolyWrapper STATIC ${Wsrc})
target_link_libraries(NTPolyWrapper NTPoly)

include(GNUInstallDirs)
install(TARGETS NTPolyWrapper
        EXPORT ntpolywrapper-export
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

################################################################################
install(EXPORT ntpolywrapper-export
        FILE NTPolyWrapperConfig.cmake
        NAMESPACE NTPoly::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NTPoly)
