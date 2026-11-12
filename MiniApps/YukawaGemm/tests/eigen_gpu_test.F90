PROGRAM EigenGPUTest
  USE ISO_C_BINDING
  USE NTMPIModule
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE EigenSolversModule, ONLY: EigenDecomposition
  USE ExponentialSolversModule, ONLY: ComputeDenseExponential
  USE DensityMatrixSolversModule, ONLY: DenseDensity
  USE SolverParametersModule
  USE TripletModule
  USE TripletListModule
  IMPLICIT NONE
  TYPE(Matrix_ps) :: h,wcpu,wgpu,qcpu,qgpu,reference,check_matrix,temp,identity,qt,pcpu,pgpu
  TYPE(TripletList_r) :: entries
  TYPE(Triplet_r) :: entry
  TYPE(TripletList_c) :: complex_entries
  TYPE(Triplet_c) :: complex_entry
  LOGICAL :: complex_problem
  TYPE(SolverParameters_t) :: params
  INTEGER :: ierr,rank,nr,pr,n,example,mode,k,j,status,slices,work_comm
  REAL(NTREAL) :: error,t0,seconds,residual,orthogonality,dc,dg,ht,ct,gt
  CHARACTER(1024) :: arg
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_CHAR,C_INT
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
  END INTERFACE
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,arg);READ(arg,*) pr
  CALL GET_COMMAND_ARGUMENT(2,arg);READ(arg,*) n
  slices=1;CALL GET_COMMAND_ARGUMENT(3,arg);IF(LEN_TRIM(arg)>0) READ(arg,*) slices
  CALL GET_COMMAND_ARGUMENT(4,arg)
  IF(TRIM(arg)=='split') THEN
    CALL MPI_Comm_split(MPI_COMM_WORLD,MOD(rank,2),rank,work_comm,ierr)
    CALL MPI_Comm_rank(work_comm,rank,ierr);CALL MPI_Comm_size(work_comm,nr,ierr)
  ELSE
    CALL MPI_Comm_dup(MPI_COMM_WORLD,work_comm,ierr)
  END IF
  CALL ConstructProcessGrid(work_comm,pr,nr/pr/slices,slices)
  CALL GET_COMMAND_ARGUMENT(5,arg);complex_problem=TRIM(arg)=='complex'
  CALL ConstructEmptyMatrix(h,n,is_complex_in=complex_problem)
  CALL ConstructSolverParameters(params,threshold_in=0.0_NTREAL,be_verbose_in=.FALSE.)
  DO example=1,8
    CALL FillInput()
    k=n
    IF(example==2) k=MIN(3,n)
    IF(example==3) k=0
    IF(example==4) k=n+2
    params%threshold=0
    IF(example==5) params%threshold=0.05_NTREAL
    DO mode=0,1
      CALL SetMode(mode)
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr);t0=MPI_Wtime()
      IF(mode==0) THEN
        IF(example==6) THEN
          CALL EigenDecomposition(h,wcpu,nvals_in=k,solver_parameters_in=params)
        ELSE
          CALL EigenDecomposition(h,wcpu,qcpu,k,params)
        END IF
      ELSE
        IF(example==6) THEN
          CALL EigenDecomposition(h,wgpu,nvals_in=k,solver_parameters_in=params)
        ELSE
          CALL EigenDecomposition(h,wgpu,qgpu,k,params)
        END IF
      END IF
      seconds=MPI_Wtime()-t0
      CALL MPI_Allreduce(MPI_IN_PLACE,seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,work_comm,ierr)
      IF(rank==0) PRINT *, 'eigen_time',example,mode,seconds
    END DO
    CALL Finite(wgpu)
    IF(example==1) THEN
      CALL MatrixTrace(h,ht);CALL MatrixTrace(wcpu,ct);CALL MatrixTrace(wgpu,gt)
      CALL MatrixMultiply(h,qcpu,temp);CALL MatrixMultiply(qcpu,wcpu,check_matrix)
      CALL IncrementMatrix(temp,check_matrix,alpha_in=-1.0_NTREAL);dc=MatrixNorm(check_matrix)
      CALL MatrixMultiply(h,qgpu,temp);CALL MatrixMultiply(qgpu,wgpu,check_matrix)
      CALL IncrementMatrix(temp,check_matrix,alpha_in=-1.0_NTREAL);dg=MatrixNorm(check_matrix)
      IF(rank==0) PRINT *, 'eigen_diagnostics input/cpu/gpu_trace cpu/gpu_residual',ht,ct,gt,dc,dg
    END IF
    CALL Difference(wcpu,wgpu,error)
    IF(.NOT.(error<1e-10_NTREAL)) CALL Fail('eigenvalues',error)
    residual=0;orthogonality=0
    IF(example/=6) THEN
      CALL Finite(qgpu)
      CALL TransposeMatrix(qcpu,qt);CALL ConjugateMatrix(qt);CALL MatrixMultiply(qcpu,qt,pcpu)
      CALL TransposeMatrix(qgpu,qt);CALL ConjugateMatrix(qt);CALL MatrixMultiply(qgpu,qt,pgpu)
      CALL Difference(pcpu,pgpu,error)
      IF(.NOT.(error<1e-9_NTREAL)) CALL Fail('projector',error)
      IF(params%threshold==0 .AND. example/=8) THEN
        CALL MatrixMultiply(h,qgpu,check_matrix);CALL MatrixMultiply(qgpu,wgpu,temp)
        CALL IncrementMatrix(temp,check_matrix,alpha_in=-1.0_NTREAL)
        residual=MatrixNorm(check_matrix)/MAX(1.0_NTREAL,MatrixNorm(h))
        IF(.NOT.(residual<1e-10_NTREAL)) CALL Fail('residual',residual)
        CALL MatrixMultiply(qt,qgpu,check_matrix)
        CALL ConstructEmptyMatrix(identity,h);CALL ConstructTripletList(entries)
        IF(rank==0) THEN
          DO j=1,MIN(k,n)
            entry%index_row=j;entry%index_column=j;entry%point_value=1
            CALL AppendToTripletList(entries,entry)
          END DO
        END IF
        IF(complex_problem) THEN
          CALL CopyTripletList(entries,complex_entries)
          CALL FillMatrixFromTripletList(identity,complex_entries);CALL DestructTripletList(complex_entries)
        ELSE
          CALL FillMatrixFromTripletList(identity,entries)
        END IF
        CALL DestructTripletList(entries)
        CALL Difference(identity,check_matrix,orthogonality)
        IF(.NOT.(orthogonality<1e-10_NTREAL)) CALL Fail('orthogonality',orthogonality)
      END IF
    END IF
    IF(rank==0) PRINT *, 'EIGEN GPU API PASS',example,error,residual,orthogonality
    CALL DestructMatrix(wcpu);CALL DestructMatrix(wgpu);CALL DestructMatrix(qcpu);CALL DestructMatrix(qgpu)
  END DO
  CALL FillInput()
  CALL ComplexFallback()
  params%threshold=1e-12_NTREAL
  CALL SetMode(0);CALL ComputeDenseExponential(h,pcpu,params)
  CALL SetMode(1);CALL ComputeDenseExponential(h,pgpu,params)
  CALL Difference(pcpu,pgpu,error);IF(.NOT.(error<1e-9_NTREAL)) CALL Fail('dense exponential',error)
  CALL ConstructEmptyMatrix(identity,h);CALL FillMatrixIdentity(identity)
  CALL SetMode(0);CALL DenseDensity(h,identity,REAL(n/2,NTREAL),pcpu,solver_parameters_in=params)
  CALL SetMode(1);CALL DenseDensity(h,identity,REAL(n/2,NTREAL),pgpu,solver_parameters_in=params)
  CALL Difference(pcpu,pgpu,error);IF(.NOT.(error<1e-9_NTREAL)) CALL Fail('dense density',error)
  IF(rank==0) PRINT *, 'DENSE CALLERS GPU PASS'
  CALL DestructMatrix(h);CALL DestructMatrix(identity);CALL DestructMatrix(reference);CALL DestructMatrix(check_matrix)
  CALL DestructMatrix(temp);CALL DestructMatrix(qt);CALL DestructMatrix(pcpu);CALL DestructMatrix(pgpu)
  CALL DestructSolverParameters(params);CALL DestructProcessGrid()
  IF(rank==0) PRINT *, 'EIGEN GPU TEARDOWN PASS'
  CALL MPI_Comm_free(work_comm,ierr);CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE SetMode(value)
    INTEGER, INTENT(IN) :: value
    CHARACTER :: text
    WRITE(text,'(I1)') value
    status=setenv('NTPOLY_GPU_EIGEN'//C_NULL_CHAR,text//C_NULL_CHAR,1_C_INT)
  END SUBROUTINE
  SUBROUTINE FillInput()
    INTEGER :: i
    CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,n
        entry%index_row=i;entry%index_column=i
        entry%point_value=1.0_NTREAL+REAL(i,NTREAL)/n+0.001_NTREAL*example
        CALL AppendToTripletList(entries,entry)
        IF(i<n) THEN
          entry%index_row=i+1;entry%index_column=i;entry%point_value=0.07_NTREAL
          IF(example==8) entry%point_value=0.19_NTREAL
          CALL AppendToTripletList(entries,entry)
          entry%index_row=i;entry%index_column=i+1;entry%point_value=0.07_NTREAL
          CALL AppendToTripletList(entries,entry)
        END IF
      END DO
    END IF
    IF(complex_problem) THEN
      CALL CopyTripletList(entries,complex_entries)
      DO i=1,complex_entries%CurrentSize
        IF(complex_entries%DATA(i)%index_row<complex_entries%DATA(i)%index_column) THEN
          complex_entries%DATA(i)%point_value=complex_entries%DATA(i)%point_value+CMPLX(0.0_NTREAL,.04_NTREAL,NTREAL)
        ELSE IF(complex_entries%DATA(i)%index_row>complex_entries%DATA(i)%index_column) THEN
          complex_entries%DATA(i)%point_value=complex_entries%DATA(i)%point_value-CMPLX(0.0_NTREAL,.04_NTREAL,NTREAL)
        END IF
        IF(example==8 .AND. complex_entries%DATA(i)%index_row==complex_entries%DATA(i)%index_column) &
          complex_entries%DATA(i)%point_value=complex_entries%DATA(i)%point_value+CMPLX(0.0_NTREAL,.2_NTREAL,NTREAL)
      END DO
      CALL FillMatrixFromTripletList(h,complex_entries);CALL DestructTripletList(complex_entries)
    ELSE
      CALL FillMatrixFromTripletList(h,entries)
    END IF
    CALL DestructTripletList(entries)
  END SUBROUTINE
  SUBROUTINE ComplexFallback()
    TYPE(Matrix_ps) :: ch,cw,gw,cq,gq,ct,p1,p2
    TYPE(TripletList_c) :: list
    TYPE(Triplet_c) :: z
    INTEGER :: i
    CALL ConstructEmptyMatrix(ch,n,process_grid_in=h%process_grid,is_complex_in=.TRUE.)
    CALL ConstructTripletList(list)
    IF(rank==0) THEN
      DO i=1,n
        z%index_row=i;z%index_column=i;z%point_value=CMPLX(1.0_NTREAL+REAL(i,NTREAL)/n,0.0_NTREAL,KIND=NTREAL)
        CALL AppendToTripletList(list,z)
        IF(i<n) THEN
          z%index_row=i;z%index_column=i+1;z%point_value=CMPLX(.03_NTREAL,.02_NTREAL,KIND=NTREAL)
          CALL AppendToTripletList(list,z)
          z%index_row=i+1;z%index_column=i;z%point_value=CONJG(z%point_value)
          CALL AppendToTripletList(list,z)
        END IF
      END DO
    END IF
    CALL FillMatrixFromTripletList(ch,list);CALL DestructTripletList(list)
    params%threshold=0
    CALL SetMode(0);CALL EigenDecomposition(ch,cw,cq,solver_parameters_in=params)
    CALL SetMode(1);CALL EigenDecomposition(ch,gw,gq,solver_parameters_in=params)
    CALL Difference(cw,gw,error);IF(.NOT.(error<1e-12_NTREAL)) CALL Fail('complex eigenvalues',error)
    CALL TransposeMatrix(cq,ct);CALL ConjugateMatrix(ct);CALL MatrixMultiply(cq,ct,p1)
    CALL TransposeMatrix(gq,ct);CALL ConjugateMatrix(ct);CALL MatrixMultiply(gq,ct,p2)
    CALL Difference(p1,p2,error);IF(.NOT.(error<1e-12_NTREAL)) CALL Fail('complex projector',error)
    CALL DestructMatrix(ch);CALL DestructMatrix(cw);CALL DestructMatrix(gw);CALL DestructMatrix(cq)
    CALL DestructMatrix(gq);CALL DestructMatrix(ct);CALL DestructMatrix(p1);CALL DestructMatrix(p2)
    IF(rank==0) PRINT *, 'COMPLEX EIGEN COMPARISON PASS'
  END SUBROUTINE
  SUBROUTINE Difference(a,b,result)
    TYPE(Matrix_ps), INTENT(IN) :: a,b
    REAL(NTREAL), INTENT(OUT) :: result
    CALL CopyMatrix(b,reference);CALL IncrementMatrix(a,reference,alpha_in=-1.0_NTREAL)
    result=MatrixNorm(reference)/MAX(1.0_NTREAL,MatrixNorm(a))
  END SUBROUTINE
  SUBROUTINE Finite(a)
    TYPE(Matrix_ps), INTENT(IN) :: a
    INTEGER :: i,j
    IF(a%is_complex) THEN
      DO j=1,SIZE(a%local_data_c,2)
        DO i=1,SIZE(a%local_data_c,1)
          IF(ANY(.NOT.IEEE_IS_FINITE(REAL(a%local_data_c(i,j)%values,NTREAL))) .OR. &
             ANY(.NOT.IEEE_IS_FINITE(AIMAG(a%local_data_c(i,j)%values)))) CALL Fail('nonfinite complex output',1.0_NTREAL)
        END DO
      END DO
      RETURN
    END IF
    DO j=1,SIZE(a%local_data_r,2)
      DO i=1,SIZE(a%local_data_r,1)
        IF(ANY(.NOT.IEEE_IS_FINITE(a%local_data_r(i,j)%values))) CALL Fail('nonfinite output',1.0_NTREAL)
      END DO
    END DO
  END SUBROUTINE
  SUBROUTINE Fail(message,value)
    CHARACTER(*), INTENT(IN) :: message
    REAL(NTREAL), INTENT(IN) :: value
    PRINT *, 'EIGEN GPU FAILURE ',message,value
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END SUBROUTINE
END PROGRAM
