PROGRAM RealWorkspaceTest
  USE ISO_C_BINDING
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE PMatrixMemoryPoolModule
  USE TripletListModule
  USE TripletModule
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: grid,other_grid
  TYPE(MatrixMemoryPool_p) :: pool,copy_pool,other_pool
  TYPE(Matrix_ps) :: a,b,cpu,gpu,difference,other,other_result,ca,cb,ccpu
  TYPE(TripletList_r) :: entries
  TYPE(Triplet_r) :: entry
  INTEGER, PARAMETER :: sizes(6)=[129,7,65,257,33,129]
  INTEGER :: ierr,rank,nr,pr,example,n,i,rep,status,expected
  INTEGER(C_INT64_T) :: allocations
  REAL(NTREAL) :: error
  CHARACTER(64) :: arg,policy
  INTERFACE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
    INTEGER(C_INT) FUNCTION cached_slots() BIND(C,NAME='gpu_ops_workspace_slot_count_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION live_slots() BIND(C,NAME='gpu_ops_live_count_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT64_T) FUNCTION allocation_count() BIND(C,NAME='gpu_workspace_allocation_count_')
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
  expected=1;IF(TRIM(policy)=='root' .AND. rank/=0) expected=0
  CALL ConstructNewProcessGrid(grid,MPI_COMM_WORLD,pr,nr/pr,1)
  DO example=1,SIZE(sizes)
    n=sizes(example)
    CALL ConstructEmptyMatrix(a,n,process_grid_in=grid)
    CALL ConstructTripletList(entries)
    IF(rank==0) THEN
      DO i=1,n
        entry%index_row=i;entry%index_column=i;entry%point_value=1.0_NTREAL+REAL(i,NTREAL)/n
        CALL AppendToTripletList(entries,entry)
        IF(i<n) THEN
          entry%index_column=i+1;entry%point_value=.0625_NTREAL
          CALL AppendToTripletList(entries,entry)
        END IF
      END DO
    END IF
    CALL FillMatrixFromTripletList(a,entries);CALL DestructTripletList(entries)
    CALL TransposeMatrix(a,b)
    DO rep=1,2
      CALL SetMode(0);CALL MatrixMultiply(a,b,cpu)
      allocations=allocation_count()
      CALL SetMode(1);CALL MatrixMultiply(a,b,gpu,memory_pool_in=pool)
      IF(rep==2 .AND. allocation_count()/=allocations) CALL Fail('warm allocation')
      CALL Compare(cpu,gpu)
      IF(cached_slots()/=4*expected .OR. live_slots()/=0) CALL Fail('pool slot ownership')
      CALL ScaleMatrix(a,.875_NTREAL)
    END DO
    IF(rank==0) PRINT *, 'REAL WORKSPACE SHAPE/FRESH/NOALLOC PASS',example,n
  END DO
  ! Assignment retains scratch safely; destroying either copy is independent.
  copy_pool=pool
  CALL DestructMatrixMemoryPool(pool)
  CALL SetMode(0);CALL MatrixMultiply(a,b,cpu)
  CALL SetMode(1);CALL MatrixMultiply(a,b,gpu,memory_pool_in=copy_pool);CALL Compare(cpu,gpu)
  IF(cached_slots()/=4*expected) CALL Fail('copied pool lifetime')
  CALL ScopedPool()
  IF(cached_slots()/=4*expected) CALL Fail('automatic pool finalization')
  ! The same pool can contain complex scratch without corrupting real buffers.
  CALL ConvertMatrixToComplex(a,ca);CALL ConvertMatrixToComplex(b,cb)
  CALL SetMode(0);CALL MatrixMultiply(ca,cb,ccpu)
  CALL SetMode(1);CALL MatrixMultiply(ca,cb,gpu,memory_pool_in=copy_pool);CALL Compare(ccpu,gpu)
  CALL DestructMatrix(ca);CALL DestructMatrix(cb);CALL DestructMatrix(ccpu)
  CALL ConstructNewProcessGrid(other_grid,MPI_COMM_WORLD,pr,nr/pr,1)
  CALL ConstructEmptyMatrix(other,17,process_grid_in=other_grid);CALL FillMatrixIdentity(other)
  CALL MatrixMultiply(other,other,other_result,memory_pool_in=other_pool)
  IF(cached_slots()/=8*expected) CALL Fail('independent pools')
  CALL DestructMatrix(other);CALL DestructMatrix(other_result);CALL DestructProcessGrid(other_grid)
  IF(cached_slots()/=8*expected) CALL Fail('pool survives old grid')
  CALL MatrixMultiply(a,b,gpu,memory_pool_in=other_pool);CALL Compare(cpu,gpu)
  CALL DestructMatrixMemoryPool(other_pool)
  IF(cached_slots()/=4*expected) CALL Fail('independent pool destruction')
  CALL DestructMatrixMemoryPool(copy_pool)
  IF(cached_slots()/=0) CALL Fail('explicit destruction')
  ! No supplied pool still gets bounded grid-owned scratch.
  CALL MatrixMultiply(a,b,gpu);CALL Compare(cpu,gpu)
  allocations=allocation_count();CALL MatrixMultiply(a,b,gpu)
  IF(allocation_count()/=allocations) CALL Fail('grid cache warm allocation')
  CALL Compare(cpu,gpu)
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL MatrixMultiply(a,b,gpu)
  IF(cached_slots()/=0) CALL Fail('cap eviction')
  CALL Compare(cpu,gpu)
  status=setenv('NTPOLY_GPU_MEMORY_MB'//C_NULL_CHAR,''//C_NULL_CHAR,1_C_INT)
  CALL MatrixMultiply(a,b,gpu,memory_pool_in=copy_pool);CALL Compare(cpu,gpu)
  CALL DestructMatrix(a);CALL DestructMatrix(b);CALL DestructMatrix(cpu)
  CALL DestructMatrix(gpu);CALL DestructMatrix(difference);CALL DestructProcessGrid(grid)
  IF(cached_slots()/=4*expected) CALL Fail('pool survives grid teardown')
  CALL DestructMatrixMemoryPool(copy_pool)
  IF(cached_slots()/=0 .OR. live_slots()/=0 .OR. context_count()/=0) CALL Fail('final teardown')
  IF(rank==0) PRINT *, 'REAL WORKSPACE POOL/COPY/FINALIZE/GRID/EVICTION PASS'
  CALL MPI_Finalize(ierr)
CONTAINS
  SUBROUTINE SetMode(mode)
    INTEGER, INTENT(IN) :: mode
    CHARACTER :: value
    WRITE(value,'(I1)') mode
    status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,value//C_NULL_CHAR,1_C_INT)
  END SUBROUTINE
  SUBROUTINE Compare(reference,result)
    TYPE(Matrix_ps), INTENT(IN) :: reference,result
    CALL CopyMatrix(result,difference);CALL IncrementMatrix(reference,difference,alpha_in=-1.0_NTREAL)
    error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(reference))
    IF(.NOT.(error<1e-11_NTREAL)) CALL Fail('numerical comparison')
    CALL DestructMatrix(difference)
  END SUBROUTINE
  SUBROUTINE ScopedPool()
    TYPE(MatrixMemoryPool_p) :: temporary,original(2),copies(2)
    CALL MatrixMultiply(a,b,gpu,memory_pool_in=temporary);CALL Compare(cpu,gpu)
    IF(cached_slots()/=8*expected) CALL Fail('scoped pool allocation')
    original=temporary
    copies=original
    CALL DestructMatrixMemoryPool(temporary)
    CALL DestructMatrixMemoryPool(original(1));CALL DestructMatrixMemoryPool(original(2))
    CALL MatrixMultiply(a,b,gpu,memory_pool_in=copies(2));CALL Compare(cpu,gpu)
    IF(cached_slots()/=8*expected) CALL Fail('array pool assignment')
  END SUBROUTINE
  SUBROUTINE Fail(message)
    CHARACTER(*), INTENT(IN) :: message
    PRINT *, 'REAL WORKSPACE FAILURE ',message,rank,example,rep,error
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END SUBROUTINE
END PROGRAM
