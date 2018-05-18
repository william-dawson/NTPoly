set(Fsrc
    ChebyshevSolversModule.f90
    DataTypesModule.f90
    DenseMatrixModule.f90
    DensityMatrixSolversModule.f90
    DistributedMatrixMemoryPoolModule.f90
    DistributedSparseMatrixModule.f90
    DistributedSparseMatrixAlgebraModule.f90
    EigenBoundsModule.f90
    EigenSolversModule.f90
    ErrorModule.f90
    ExponentialSolversModule.f90
    FixedSolversModule.f90
    GemmTasksModule.f90
    GeometryOptimizationModule.f90
    HermiteSolversModule.f90
    InverseSolversModule.f90
    JacobiEigenSolverModule.f90
    IterativeSolversModule.f90
    LinearSolversModule.f90
    LoadBalancerModule.f90
    LoggingModule.f90
    MatrixBroadcastModule.f90
    MatrixGatherModule.f90
    MatrixMarketModule.f90
    MatrixMemoryPoolModule.f90
    MatrixSendRecvModule.f90
    MinimizerSolversModule.f90
    ParameterConverterModule.f90
    PermutationModule.f90
    PolynomialSolversModule.f90
    ProcessGridModule.f90
    RootSolversModule.f90
    SignSolversModule.f90
    SparseMatrixAlgebraModule.f90
    SparseMatrixIteratorModule.f90
    SparseMatrixModule.f90
    SparseVectorModule.f90
    SquareRootSolversModule.f90
    TimerModule.f90
    TrigonometrySolversModule.f90
    TripletListModule.f90
    TripletModule.f90
)

add_library(NTPoly STATIC ${Fsrc})
install(TARGETS NTPoly
        ARCHIVE DESTINATION lib)
