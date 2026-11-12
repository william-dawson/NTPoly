! CPU/GPU comparisons through the public algebra and iterative solver APIs.
PROGRAM AlgebraGPUTest
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE ISO_C_BINDING
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE SMatrixModule, ONLY: Matrix_lsr, DestructLocalMatrix=>DestructMatrix
  USE PSMatrixAlgebraModule
  USE SolverParametersModule
  USE InverseSolversModule
  USE SquareRootSolversModule
  USE RootSolversModule
  USE SignSolversModule
  USE ExponentialSolversModule
  USE TrigonometrySolversModule
  USE LinearSolversModule, ONLY: CGSolver, CholeskyDecomposition
  USE EigenSolversModule, ONLY: EigenDecomposition, EstimateGap
  USE SingularValueSolversModule
  USE GeometryOptimizationModule
  USE DensityMatrixSolversModule
  USE FermiOperatorModule
  USE PolynomialSolversModule, ONLY: Polynomial_t, PC=>ConstructPolynomial, PD=>DestructPolynomial, &
    PS=>SetCoefficient, PolynomialCompute=>Compute, PolynomialFactorized=>FactorizedCompute
  USE ChebyshevSolversModule, ONLY: ChebyshevPolynomial_t, CC=>ConstructPolynomial, CD=>DestructPolynomial, &
    CS=>SetCoefficient, ChebyshevCompute=>Compute, ChebyshevFactorized=>FactorizedCompute
  USE HermiteSolversModule, ONLY: HermitePolynomial_t, HC=>ConstructPolynomial, HD=>DestructPolynomial, &
    HS=>SetCoefficient, HermiteCompute=>Compute
  USE TripletModule, ONLY: Triplet_r, Triplet_c
  USE TripletListModule
  IMPLICIT NONE
  TYPE(Matrix_ps) :: a,b,h,identity,cpu,gpu,difference,aux,aux2
  TYPE(SolverParameters_t) :: params
  TYPE(Polynomial_t) :: poly
  TYPE(ChebyshevPolynomial_t) :: cheby
  TYPE(HermitePolynomial_t) :: hermite
  INTEGER :: ierr,rank,nr,pr,case_id,j,k,mode,status,slices
  REAL(NTREAL) :: alpha,beta,threshold,error,ec,eg,tc,tg
  CHARACTER(32) :: arg
  LOGICAL :: complex_problem
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
    INTEGER(C_INT) FUNCTION live_count() BIND(C,NAME='gpu_ops_live_count_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION context_count() BIND(C,NAME='solver_nccl_context_count_')
      IMPORT C_INT
    END FUNCTION
  END INTERFACE
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,arg); READ(arg,*) pr
  slices=1
  CALL GET_COMMAND_ARGUMENT(2,arg)
  IF(LEN_TRIM(arg)>0) READ(arg,*) slices
  CALL GET_COMMAND_ARGUMENT(3,arg);complex_problem=TRIM(arg)=='complex'
  CALL ConstructProcessGrid(MPI_COMM_WORLD,pr,nr/pr/slices,slices)
  status=setenv('NTPOLY_GPU_TRS2'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL ConstructEmptyMatrix(a,8,is_complex_in=complex_problem); CALL ConstructEmptyMatrix(b,8,is_complex_in=complex_problem)
  CALL ConstructEmptyMatrix(h,8,is_complex_in=complex_problem); CALL ConstructEmptyMatrix(identity,8,is_complex_in=complex_problem)
  CALL FillMatrixIdentity(identity)
  tc=0; tg=0
  CALL FillInputs()
  CALL ConstructSolverParameters(params,converge_diff_in=1e-10_NTREAL,threshold_in=1e-12_NTREAL, &
    max_iterations_in=100,be_verbose_in=.FALSE.)
  CALL PC(poly,8); CALL CC(cheby,8); CALL HC(hermite,8)
  ! The Fortran coefficient API uses positions 1..degree (constant first).
  DO j=1,8
    CALL PS(poly,j,0.01_NTREAL/(j+1)); CALL CS(cheby,j,0.01_NTREAL/(j+1))
    CALL HS(hermite,j,0.01_NTREAL/(j+1))
  END DO
  DO case_id=1,18
    alpha=1; beta=0; threshold=0
    SELECT CASE(case_id)
    CASE(2); alpha=-1
    CASE(3); alpha=2; threshold=0.1_NTREAL
    CASE(4); alpha=0.5_NTREAL; threshold=0.1_NTREAL
    CASE(5); alpha=0
    CASE(6); beta=1
    CASE(7); beta=-0.5_NTREAL; threshold=0.1_NTREAL
    CASE(8); alpha=-2; beta=2; threshold=0.1_NTREAL
    CASE(9); beta=TINY(beta)/2
    CASE(10); alpha=2; threshold=1.0_NTREAL
    CASE(11); CALL ConstructEmptyMatrix(b,8,is_complex_in=complex_problem)
    CASE(12); CALL FillInputs()
    CASE(13); CALL ConstructEmptyMatrix(a,8,is_complex_in=complex_problem)
    CASE(14); CALL FillInputs(); CALL ScaleMatrix(b,-0.7_NTREAL)
    CASE(15); alpha=0; beta=1
    CASE(16); threshold=1e3_NTREAL; beta=1
    CASE(17); alpha=-0.5_NTREAL; beta=-1; threshold=0.02_NTREAL
    CASE(18); alpha=1e-3_NTREAL; beta=2; threshold=1e-3_NTREAL
    END SELECT
    DO mode=0,1
      CALL SetMode(mode)
      IF(mode==0) THEN
        CALL CopyMatrix(b,cpu)
        CALL MatrixMultiply(a,b,cpu,alpha,beta,threshold)
      ELSE
        CALL CopyMatrix(b,gpu)
        CALL MatrixMultiply(a,b,gpu,alpha,beta,threshold)
      END IF
    END DO
    CALL Compare('ALGEBRA',case_id)
  END DO
  CALL FillInputs()
  CALL ComplexFallbackTest()
  CALL ComplexEdgeTest()
  CALL PaddedGatherTest()
  DO case_id=1,49
    IF(complex_problem .AND. case_id==42) THEN
      IF(rank==0) PRINT *, 'SOLVER SKIP real-only Cholesky',case_id
      CYCLE
    END IF
    DO mode=0,1
      CALL SetMode(mode)
      CALL DestructMatrix(aux); CALL DestructMatrix(aux2)
      IF(mode==0) THEN
        CALL Solve(cpu)
      ELSE
        CALL Solve(gpu)
      END IF
    END DO
    CALL Compare('SOLVER',case_id)
  END DO
  CALL DestructMatrix(a); CALL DestructMatrix(b); CALL DestructMatrix(h)
  CALL DestructMatrix(identity); CALL DestructMatrix(aux); CALL DestructMatrix(aux2)
  CALL PD(poly); CALL CD(cheby); CALL HD(hermite); CALL DestructSolverParameters(params)
  CALL DestructProcessGrid()
  IF(context_count()/=0) CALL MPI_Abort(MPI_COMM_WORLD,22,ierr)
  IF(rank==0) PRINT *, 'ALGEBRA SOLVER COVERAGE COMPLETE'
  CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE SetMode(m)
    INTEGER, INTENT(IN) :: m
    CHARACTER :: value
    WRITE(value,'(I1)') m
    status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,value//C_NULL_CHAR,1_C_INT)
    status=setenv('NTPOLY_GPU_EIGEN'//C_NULL_CHAR,value//C_NULL_CHAR,1_C_INT)
  END SUBROUTINE
  SUBROUTINE FillInputs()
    TYPE(TripletList_r) :: entries
    TYPE(Triplet_r) :: entry
    INTEGER :: i,l
    CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,8
        entry%index_row=i; entry%index_column=i; entry%point_value=1.0_NTREAL+0.1_NTREAL*i
        CALL AppendToTripletList(entries,entry)
        IF(i<8) THEN
          entry%index_column=i+1; entry%point_value=0.0625_NTREAL
          CALL AppendToTripletList(entries,entry)
          entry%index_row=i+1; entry%index_column=i
          CALL AppendToTripletList(entries,entry)
        END IF
      END DO
    END IF
    CALL FillTyped(a,entries,.TRUE.)
    CALL DestructTripletList(entries); CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,8
        DO l=1,8
          IF(MOD(i+2*l,3)/=0) CYCLE
          entry%index_row=i; entry%index_column=l; entry%point_value=(i-2*l)*0.03125_NTREAL
          CALL AppendToTripletList(entries,entry)
        END DO
      END DO
    END IF
    CALL FillTyped(b,entries,.FALSE.)
    CALL DestructTripletList(entries); CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,8
        entry%index_row=i; entry%index_column=i
        entry%point_value=SIGN(1.0_NTREAL,REAL(i,NTREAL)-4.5_NTREAL)
        CALL AppendToTripletList(entries,entry)
        IF(i<8) THEN
          entry%index_column=i+1; entry%point_value=0.04_NTREAL
          CALL AppendToTripletList(entries,entry)
          entry%index_row=i+1; entry%index_column=i
          CALL AppendToTripletList(entries,entry)
        END IF
      END DO
    END IF
    CALL FillTyped(h,entries,.TRUE.); CALL DestructTripletList(entries)
  END SUBROUTINE
  SUBROUTINE FillTyped(matrix,entries,hermitian)
    TYPE(Matrix_ps), INTENT(INOUT) :: matrix
    TYPE(TripletList_r), INTENT(IN) :: entries
    LOGICAL, INTENT(IN) :: hermitian
    TYPE(TripletList_c) :: list
    INTEGER :: i,r,c
    REAL(NTREAL) :: imag
    IF(.NOT.complex_problem) THEN
      CALL FillMatrixFromTripletList(matrix,entries)
      RETURN
    END IF
    CALL CopyTripletList(entries,list)
    DO i=1,list%CurrentSize
      r=list%DATA(i)%index_row;c=list%DATA(i)%index_column
      imag=.015625_NTREAL*(r-c)
      IF(.NOT.hermitian) imag=.0078125_NTREAL*(r+2*c)
      list%DATA(i)%point_value=list%DATA(i)%point_value+CMPLX(0.0_NTREAL,imag,NTREAL)
    END DO
    CALL FillMatrixFromTripletList(matrix,list);CALL DestructTripletList(list)
  END SUBROUTINE
  SUBROUTINE Compare(kind,id)
    CHARACTER(*), INTENT(IN) :: kind
    INTEGER, INTENT(IN) :: id
    CALL CheckFinite(cpu);CALL CheckFinite(gpu)
    CALL CopyMatrix(gpu,difference)
    CALL IncrementMatrix(cpu,difference,alpha_in=-1.0_NTREAL)
    error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(cpu))
    IF(.NOT.(error<1e-8_NTREAL)) THEN
      PRINT *, kind,' PARITY FAILURE',id,error
      CALL MPI_Abort(MPI_COMM_WORLD,20,ierr)
    END IF
    IF(live_count()/=0) CALL MPI_Abort(MPI_COMM_WORLD,21,ierr)
    IF(rank==0) PRINT *, kind,' PARITY PASS',id,error,tc,tg
    CALL DestructMatrix(cpu); CALL DestructMatrix(gpu); CALL DestructMatrix(difference)
  END SUBROUTINE
  SUBROUTINE CheckFinite(matrix)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    INTEGER :: i,j
    IF(matrix%is_complex) THEN
      DO j=1,SIZE(matrix%local_data_c,2)
        DO i=1,SIZE(matrix%local_data_c,1)
          IF(ANY(.NOT.IEEE_IS_FINITE(REAL(matrix%local_data_c(i,j)%values,NTREAL))) .OR. &
             ANY(.NOT.IEEE_IS_FINITE(AIMAG(matrix%local_data_c(i,j)%values)))) CALL MPI_Abort(MPI_COMM_WORLD,27,ierr)
        END DO
      END DO
    ELSE
      DO j=1,SIZE(matrix%local_data_r,2)
        DO i=1,SIZE(matrix%local_data_r,1)
          IF(ANY(.NOT.IEEE_IS_FINITE(matrix%local_data_r(i,j)%values))) CALL MPI_Abort(MPI_COMM_WORLD,26,ierr)
        END DO
      END DO
    END IF
  END SUBROUTINE
  SUBROUTINE PaddedGatherTest()
    TYPE(Matrix_ps) :: padded
    TYPE(Matrix_lsr) :: gathered
    INTEGER, ALLOCATABLE :: permutation(:)
    INTEGER :: i,n
    CALL ConstructEmptyMatrix(padded,5)
    n=padded%logical_matrix_dimension
    ALLOCATE(permutation(n))
    DO i=1,n
      permutation(i)=i
    END DO
    CALL FillMatrixPermutation(padded,permutation)
    CALL GatherMatrixToProcess(padded,gathered,0,include_padding_in=.TRUE.)
    IF(padded%process_grid%within_slice_rank==0) THEN
      IF(gathered%rows/=n .OR. gathered%columns/=n) CALL MPI_Abort(MPI_COMM_WORLD,29,ierr)
      IF(SIZE(gathered%values)/=n) CALL MPI_Abort(MPI_COMM_WORLD,30,ierr)
      IF(ANY(gathered%values/=1.0_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,31,ierr)
    END IF
    CALL DestructLocalMatrix(gathered); CALL DestructMatrix(padded)
    IF(rank==0) PRINT *, 'PADDED LOGICAL GATHER PASS'
  END SUBROUTINE
  SUBROUTINE ComplexFallbackTest()
    TYPE(Matrix_ps) :: ca,cb,co,cg
    IF(a%is_complex) THEN
      CALL CopyMatrix(a,ca)
    ELSE
      CALL ConvertMatrixToComplex(a,ca)
    END IF
    IF(b%is_complex) THEN
      CALL CopyMatrix(b,cb)
    ELSE
      CALL ConvertMatrixToComplex(b,cb)
    END IF
    CALL ScaleMatrix(ca,CMPLX(1.0_NTREAL,0.3_NTREAL,KIND=NTREAL))
    CALL SetMode(0); CALL MatrixMultiply(ca,cb,co,threshold_in=1e-12_NTREAL)
    CALL SetMode(1); CALL MatrixMultiply(ca,cb,cg,threshold_in=1e-12_NTREAL)
    CALL IncrementMatrix(co,cg,alpha_in=-1.0_NTREAL)
    IF(.NOT.(MatrixNorm(cg)<1e-12_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,28,ierr)
    CALL DestructMatrix(ca); CALL DestructMatrix(cb); CALL DestructMatrix(co); CALL DestructMatrix(cg)
    IF(rank==0) PRINT *, 'COMPLEX PRODUCT COMPARISON PASS'
  END SUBROUTINE
  SUBROUTINE ComplexEdgeTest()
    TYPE(Matrix_ps) :: ca,cb,real_a,real_b,left,right,initial,co,cg,tmp
    TYPE(TripletList_c) :: list
    TYPE(Triplet_c) :: z
    INTEGER :: i,test
    ! Odd dimension exercises padding. A complex magnitude can pass pruning
    ! even when both its real and imaginary components are below the cutoff.
    CALL ConstructEmptyMatrix(ca,5,is_complex_in=.TRUE.)
    CALL ConstructEmptyMatrix(cb,5,is_complex_in=.TRUE.)
    CALL ConstructTripletList(list)
    IF(rank==0) THEN
      DO i=1,5
        z%index_row=i;z%index_column=i;z%point_value=CMPLX(.06_NTREAL,.08_NTREAL,NTREAL)
        CALL AppendToTripletList(list,z)
      END DO
    END IF
    CALL FillMatrixFromTripletList(ca,list);CALL DestructTripletList(list)
    CALL FillMatrixIdentity(cb)
    CALL ConvertMatrixToReal(ca,real_a);CALL ConvertMatrixToReal(cb,real_b)
    DO test=1,6
      CALL CopyMatrix(ca,left);CALL CopyMatrix(cb,right);CALL CopyMatrix(ca,initial)
      IF(test==2) CALL CopyMatrix(real_a,left)
      IF(test==3) CALL CopyMatrix(real_b,right)
      IF(test==4) CALL CopyMatrix(real_b,initial)
      alpha=1;beta=0;threshold=.09_NTREAL
      IF(test==2) threshold=.01_NTREAL
      IF(test==4) beta=-.5_NTREAL
      IF(test==5) THEN
        beta=-1;threshold=0
      END IF
      IF(test==6) THEN
        alpha=2;beta=.25_NTREAL;threshold=.15_NTREAL
      END IF
      CALL CopyMatrix(initial,co);CALL CopyMatrix(initial,cg)
      CALL SetMode(0);CALL MatrixMultiply(left,right,co,alpha,beta,threshold)
      CALL SetMode(1);CALL MatrixMultiply(left,right,cg,alpha,beta,threshold)
      CALL CopyMatrix(cg,tmp);CALL IncrementMatrix(co,tmp,alpha_in=-1.0_NTREAL)
      IF(.NOT.(MatrixNorm(tmp)<1e-12_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,32,ierr)
      IF(test==1 .AND. .NOT.(MatrixNorm(cg)>.09_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,33,ierr)
      IF(test==5 .AND. .NOT.(MatrixNorm(cg)<1e-12_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,34,ierr)
      ! Same input storage, changed values on the next call must be observed.
      CALL ScaleMatrix(left,CMPLX(.5_NTREAL,-.25_NTREAL,NTREAL))
      CALL SetMode(0);CALL MatrixMultiply(left,right,co,threshold_in=0.0_NTREAL)
      CALL SetMode(1);CALL MatrixMultiply(left,right,cg,threshold_in=0.0_NTREAL)
      CALL CopyMatrix(cg,tmp);CALL IncrementMatrix(co,tmp,alpha_in=-1.0_NTREAL)
      IF(.NOT.(MatrixNorm(tmp)<1e-12_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,35,ierr)
    END DO
    CALL DestructMatrix(ca);CALL DestructMatrix(cb);CALL DestructMatrix(real_a);CALL DestructMatrix(real_b)
    CALL DestructMatrix(left);CALL DestructMatrix(right);CALL DestructMatrix(initial)
    CALL DestructMatrix(co);CALL DestructMatrix(cg);CALL DestructMatrix(tmp)
    IF(rank==0) PRINT *, 'COMPLEX MIXED MAGNITUDE PADDING CANCELLATION FRESH INPUT PASS'
  END SUBROUTINE
  SUBROUTINE Solve(output)
    TYPE(Matrix_ps), INTENT(INOUT) :: output
    REAL(NTREAL) :: started,elapsed
    started=MPI_Wtime()
    SELECT CASE(case_id)
    CASE(1); CALL Invert(a,output,params)
    CASE(2); CALL PseudoInverse(a,output,params)
    CASE(3); CALL SquareRoot(a,output,params)
    CASE(4); CALL InverseSquareRoot(a,output,params)
    CASE(5); CALL SquareRoot(a,output,params,order_in=4)
    CASE(6); CALL InverseSquareRoot(a,output,params,order_in=4)
    CASE(7); CALL ComputeRoot(a,output,3,params)
    CASE(8); CALL ComputeInverseRoot(a,output,3,params)
    CASE(9); CALL SignFunction(h,output,params)
    CASE(10); CALL PolarDecomposition(a,output,aux,params)
    CASE(11); CALL ComputeExponential(a,output,params)
    CASE(12); CALL ComputeExponentialPade(a,output,params)
    CASE(13); CALL ComputeExponentialTaylor(a,output,params)
    CASE(14); CALL ComputeLogarithm(a,output,params)
    CASE(15); CALL ComputeLogarithmTaylor(a,output,params)
    CASE(16); CALL Sine(a,output,params)
    CASE(17); CALL Cosine(a,output,params)
    CASE(18); CALL ScaleSquareTrigonometryTaylor(a,output,params)
    CASE(19); CALL CGSolver(a,output,identity,params)
    CASE(20); CALL PolynomialCompute(a,output,poly,params)
    CASE(21); CALL PolynomialFactorized(a,output,poly,params)
    CASE(22); CALL ChebyshevCompute(a,output,cheby,params)
    CASE(23); CALL ChebyshevFactorized(a,output,cheby,params)
    CASE(24); CALL HermiteCompute(a,output,hermite,params)
    CASE(25); CALL PM(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
    CASE(26); CALL TRS2(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
    CASE(27); CALL TRS4(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
    CASE(28); CALL HPCP(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
    CASE(29); CALL ScaleAndFold(h,identity,4.0_NTREAL,output,-0.8_NTREAL,0.8_NTREAL,solver_parameters_in=params)
    CASE(30); CALL WOM_GC(h,identity,output,0.0_NTREAL,2.0_NTREAL,solver_parameters_in=params)
    CASE(31); CALL WOM_C(h,identity,output,4.0_NTREAL,2.0_NTREAL,solver_parameters_in=params)
    CASE(32); CALL DenseInvert(a,output,params)
    CASE(33); CALL DenseSquareRoot(a,output,params)
    CASE(34); CALL DenseInverseSquareRoot(a,output,params)
    CASE(35); CALL DenseSignFunction(h,output,params)
    CASE(36); CALL DenseSine(a,output,params)
    CASE(37); CALL DenseCosine(a,output,params)
    CASE(38); CALL ComputeDenseExponential(a,output,params)
    CASE(39); CALL ComputeDenseLogarithm(a,output,params)
    CASE(40); CALL DenseDensity(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
    CASE(41); CALL ComputeDenseFOE(h,identity,4.0_NTREAL,output,2.0_NTREAL,solver_parameters_in=params)
    CASE(42); CALL CholeskyDecomposition(a,output,params)
    CASE(43); CALL SingularValueDecomposition(a,aux,aux2,output,params)
    CASE(44); CALL EigenDecomposition(a,output,aux,solver_parameters_in=params)
    CASE(45); CALL EnergyDensityMatrix(h,a,output,params%threshold)
    CASE(46); CALL McWeenyStep(a,output,identity,params%threshold)
    CASE(47)
      CALL TRS2(h,identity,4.0_NTREAL,output,solver_parameters_in=params)
      CALL EstimateGap(h,output,0.0_NTREAL,elapsed,params)
      IF(mode==0) ec=elapsed
      IF(mode==1) eg=elapsed
      IF(mode==1 .AND. ABS(ec-eg)>1e-8_NTREAL) CALL MPI_Abort(MPI_COMM_WORLD,25,ierr)
    CASE(48)
      CALL TRS2(h,identity,4.0_NTREAL,aux,solver_parameters_in=params)
      CALL PurificationExtrapolate(aux,identity,4.0_NTREAL,output,params)
    CASE(49)
      CALL TRS2(h,identity,4.0_NTREAL,aux,solver_parameters_in=params)
      CALL LowdinExtrapolate(aux,identity,a,output,params)
    END SELECT
    elapsed=MPI_Wtime()-started
    IF(mode==0) THEN
      tc=elapsed
    ELSE
      tg=elapsed
    END IF
  END SUBROUTINE
END PROGRAM
