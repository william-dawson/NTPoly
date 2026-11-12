PROGRAM ComplexWorkspaceTest
  USE ISO_C_BINDING
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE TripletListModule
  USE TripletModule
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: grid,other_grid
  TYPE(Matrix_ps) :: a,b,cpu,gpu,difference,other,other_result,real_other
  TYPE(TripletList_c) :: entries
  TYPE(Triplet_c) :: entry
  INTEGER :: ierr,rank,nr,pr,example,n,i,rep,status,expected
  INTEGER, PARAMETER :: sizes(6)=[129,7,65,257,33,129]
  INTEGER(C_INT64_T) :: allocations
  REAL(NTREAL) :: error
  CHARACTER(64) :: arg,policy
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
    INTEGER(C_INT) FUNCTION workspace_count() BIND(C,NAME='ntpoly_complex_workspace_count_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT64_T) FUNCTION allocation_count() BIND(C,NAME='ntpoly_complex_allocation_count_')
      IMPORT C_INT64_T
    END FUNCTION
    INTEGER(C_INT) FUNCTION context_count() BIND(C,NAME='solver_nccl_context_count_')
      IMPORT C_INT
    END FUNCTION
  END INTERFACE
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,arg);READ(arg,*) pr
  CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_RANK_POLICY',policy)
  expected=1
  IF(TRIM(policy)=='root' .AND. rank/=0) expected=0
  CALL ConstructNewProcessGrid(grid,MPI_COMM_WORLD,pr,nr/pr,1)
  DO example=1,SIZE(sizes)
    n=sizes(example)
    CALL ConstructEmptyMatrix(a,n,process_grid_in=grid,is_complex_in=.TRUE.)
    CALL ConstructEmptyMatrix(b,a)
    CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,n
        entry%index_row=i;entry%index_column=i
        entry%point_value=CMPLX(1.0_NTREAL+REAL(i,NTREAL)/n,.125_NTREAL,NTREAL)
        CALL AppendToTripletList(entries,entry)
        IF(i<n) THEN
          entry%index_row=i;entry%index_column=i+1;entry%point_value=CMPLX(.0625_NTREAL,-.03125_NTREAL,NTREAL)
          CALL AppendToTripletList(entries,entry)
        END IF
      END DO
    END IF
    CALL FillMatrixFromTripletList(a,entries);CALL DestructTripletList(entries)
    CALL CopyMatrix(a,b);CALL TransposeMatrix(b,gpu);CALL CopyMatrix(gpu,b);CALL DestructMatrix(gpu)
    DO rep=1,2
      CALL SetMode(0);CALL MatrixMultiply(a,b,cpu,threshold_in=0.0_NTREAL)
      allocations=allocation_count()
      CALL SetMode(1);CALL MatrixMultiply(a,b,gpu,threshold_in=0.0_NTREAL)
      IF(rep==2 .AND. allocation_count()/=allocations) CALL Fail('unexpected warm device allocation')
      CALL Compare()
      IF(workspace_count()/=expected) CALL Fail('workspace owner count')
      CALL ScaleMatrix(a,CMPLX(.875_NTREAL,.0625_NTREAL,NTREAL))
    END DO
    IF(rank==0) PRINT *, 'COMPLEX WORKSPACE SHAPE/FRESH/NOALLOC PASS',example,n
  END DO
  ! A second grid must own independent scratch and release just its resources.
  CALL ConstructNewProcessGrid(other_grid,MPI_COMM_WORLD,pr,nr/pr,1)
  CALL ConstructEmptyMatrix(other,17,process_grid_in=other_grid,is_complex_in=.TRUE.)
  CALL FillMatrixIdentity(other);CALL MatrixMultiply(other,other,other_result)
  IF(workspace_count()/=2*expected) CALL Fail('independent grid workspace')
  CALL ConvertMatrixToReal(other,real_other)
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL MatrixMultiply(real_other,real_other,other_result)
  IF(workspace_count()/=0 .OR. context_count()/=2*expected) CALL Fail('cross-type scratch reclamation')
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,''//C_NULL_CHAR,1_C_INT)
  CALL SetMode(0);CALL MatrixMultiply(a,b,cpu)
  CALL SetMode(1);CALL MatrixMultiply(a,b,gpu);CALL Compare()
  CALL MatrixMultiply(other,other,other_result)
  IF(workspace_count()/=2*expected) CALL Fail('independent workspace recreation')
  IF(rank==0) PRINT *, 'COMPLEX WORKSPACE CROSS-TYPE RECLAMATION PASS'
  CALL DestructMatrix(real_other)
  CALL DestructMatrix(other);CALL DestructMatrix(other_result);CALL DestructProcessGrid(other_grid)
  IF(workspace_count()/=expected) CALL Fail('independent grid teardown')
  ! Changing the cap must evict retained scratch and take the valid CPU path.
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL MatrixMultiply(a,b,gpu)
  IF(workspace_count()/=0) CALL Fail('memory-cap eviction')
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,''//C_NULL_CHAR,1_C_INT)
  CALL SetMode(0);CALL MatrixMultiply(a,b,cpu)
  CALL SetMode(1);CALL MatrixMultiply(a,b,gpu);CALL Compare()
  IF(workspace_count()/=expected) CALL Fail('workspace recreation')
  CALL DestructMatrix(a);CALL DestructMatrix(b);CALL DestructMatrix(cpu)
  CALL DestructMatrix(gpu);CALL DestructMatrix(difference);CALL DestructProcessGrid(grid)
  IF(workspace_count()/=0 .OR. context_count()/=0) CALL Fail('workspace/context teardown')
  IF(rank==0) PRINT *, 'COMPLEX WORKSPACE EVICTION/ISOLATION/TEARDOWN PASS'
  CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE SetMode(mode)
    INTEGER, INTENT(IN) :: mode
    CHARACTER :: value
    WRITE(value,'(I1)') mode
    status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,value//C_NULL_CHAR,1_C_INT)
  END SUBROUTINE
  SUBROUTINE Compare()
    CALL CopyMatrix(gpu,difference);CALL IncrementMatrix(cpu,difference,alpha_in=-1.0_NTREAL)
    error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(cpu))
    IF(.NOT.(error<1e-11_NTREAL)) CALL Fail('numerical comparison')
    CALL DestructMatrix(difference)
  END SUBROUTINE
  SUBROUTINE Fail(message)
    CHARACTER(*), INTENT(IN) :: message
    PRINT *, 'COMPLEX WORKSPACE FAILURE ',message,rank,example,rep,error
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END SUBROUTINE
END PROGRAM
