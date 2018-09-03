################################################################################
set(Csrc
    ChebyshevSolvers.cc
    DensityMatrixSolvers.cc
    EigenBounds.cc
    EigenSolvers.cc
    ExponentialSolvers.cc
    GeometryOptimization.cc
    HermiteSolvers.cc
    InverseSolvers.cc
    LinearSolvers.cc
    LoadBalancer.cc
    MatrixMemoryPool.cc
    Permutation.cc
    Polynomial.cc
    PMatrixMemoryPool.cc
    ProcessGrid.cc
    PSMatrix.cc
    RootSolvers.cc
    SignSolvers.cc
    SMatrix.cc
    SolverBase.cc
    SolverParameters.cc
    SquareRootSolvers.cc
    TrigonometrySolvers.cc
    TripletList.cc
)

set(Chead
    ChebyshevSolvers.h
    DensityMatrixSolvers.h
    EigenBounds.h
    EigenSolvers.h
    ExponentialSolvers.h
    GeometryOptimization.h
    HermiteSolvers.h
    InverseSolvers.h
    LinearSolvers.h
    LoadBalancer.h
    MatrixMemoryPool.h
    Permutation.h
    PMatrixMemoryPool.h
    Polynomial.h
    ProcessGrid.h
    PSMatrix.h
    RootSolvers.h
    SignSolvers.h
    SMatrix.h
    SolverBase.h
    SolverParameters.h
    SquareRootSolvers.h
    TrigonometrySolvers.h
    TripletList.h
    Triplet.h
)

################################################################################
add_library(NTPolyCPP STATIC ${Csrc} ${Chead})
target_link_libraries(NTPolyCPP NTPolyWrapper)

set_target_properties(NTPolyCPP PROPERTIES PUBLIC_HEADER "${Chead}")
target_include_directories(NTPolyCPP PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Source/C>
                           $<INSTALL_INTERFACE:include>)
include(GNUInstallDirs)
install(TARGETS NTPolyCPP
        EXPORT ntpolycpp-export
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

################################################################################
# For Exporting
install(EXPORT ntpolycpp-export
        FILE NTPolyCPPConfig.cmake
        NAMESPACE NTPoly::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NTPoly)
