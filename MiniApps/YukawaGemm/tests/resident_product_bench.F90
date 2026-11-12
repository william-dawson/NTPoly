! Repeated production GPU products; input import and final validation are untimed.
PROGRAM ResidentProductBench
  USE ISO_C_BINDING
  USE OMP_LIB
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL,MPINTLONG
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule
  USE PSMatrixGPUBackendModule
  USE TripletListModule
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: grid
  TYPE(Matrix_ps) :: input,a,reference,result,difference
  TYPE(GPUMatrix_t) :: ga,gresult
  TYPE(TripletList_r) :: real_entries
  TYPE(TripletList_c) :: entries
  CHARACTER(1024) :: filename,arg
  INTEGER :: ierr,rank,nr,pr,n,kind,reps,rep,i,status,provided
  INTEGER(C_INT) :: rows,columns,dimension,row_offset,column_offset,nnz,count,profiling
  INTEGER(C_INT), ALLOCATABLE :: ir(:),ic(:),outrows(:),outcols(:)
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: iv(:),outvalues(:)
  REAL(C_DOUBLE), ALLOCATABLE :: seconds(:),phases(:,:),local_seconds(:)
  REAL(NTREAL) :: threshold,t0,error,theta
  INTEGER(KIND=8) :: input_nnz,output_nnz
  LOGICAL :: allowed,root_mode
  TYPE(C_PTR) :: output
  INTERFACE
    INTEGER(C_INT) FUNCTION resident_complex(rows,columns,n,row_offset,column_offset,nnz,ir,ic,iv, &
        threshold,reps,profile,seconds,phases,output,count) BIND(C,NAME='ntpoly_complex_resident_benchmark_')
      IMPORT C_INT,C_DOUBLE,C_DOUBLE_COMPLEX,C_PTR
      INTEGER(C_INT) :: rows,columns,n,row_offset,column_offset,nnz,reps,profile,count
      INTEGER(C_INT) :: ir(*),ic(*)
      COMPLEX(C_DOUBLE_COMPLEX) :: iv(*)
      REAL(C_DOUBLE) :: threshold,seconds(*),phases(*)
      TYPE(C_PTR) :: output
    END FUNCTION
    SUBROUTINE complex_export(output,rows,columns,values) BIND(C,NAME='ntpoly_complex_product_export_')
      IMPORT C_INT,C_DOUBLE_COMPLEX,C_PTR
      TYPE(C_PTR) :: output
      INTEGER(C_INT) :: rows(*),columns(*)
      COMPLEX(C_DOUBLE_COMPLEX) :: values(*)
    END SUBROUTINE
    SUBROUTINE profile_begin() BIND(C,NAME='gpu_profile_begin_')
    END SUBROUTINE
    SUBROUTINE profile_end() BIND(C,NAME='gpu_profile_end_')
    END SUBROUTINE
    INTEGER(C_INT) FUNCTION setenv(name,value,overwrite) BIND(C)
      IMPORT C_INT,C_CHAR
      CHARACTER(C_CHAR) :: name(*),value(*)
      INTEGER(C_INT), VALUE :: overwrite
    END FUNCTION
  END INTERFACE
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED,provided,ierr)
  IF(provided<MPI_THREAD_SERIALIZED) CALL MPI_Abort(MPI_COMM_WORLD,9,ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  CALL GET_COMMAND_ARGUMENT(1,filename)
  CALL GET_COMMAND_ARGUMENT(2,arg);READ(arg,*) pr
  CALL GET_COMMAND_ARGUMENT(3,arg);READ(arg,*) n
  CALL GET_COMMAND_ARGUMENT(4,arg);READ(arg,*) kind
  CALL GET_COMMAND_ARGUMENT(5,arg);READ(arg,*) reps
  CALL GET_COMMAND_ARGUMENT(6,arg);profiling=0;IF(TRIM(arg)=='profile') profiling=1
  IF(kind<0 .OR. kind>2 .OR. reps<2) CALL MPI_Abort(MPI_COMM_WORLD,8,ierr)
  CALL OMP_SET_NUM_THREADS(1)
  status=setenv('NTPOLY_GPU_ALGEBRA'//C_NULL_CHAR,'0'//C_NULL_CHAR,1_C_INT)
  CALL ConstructNewProcessGrid(grid,MPI_COMM_WORLD,pr,nr/pr,1)
  CALL ConstructMatrixFromBinary(input,TRIM(filename),process_grid_in=grid)
  IF(n>0 .AND. n/=input%actual_matrix_dimension) CALL ResizeMatrix(input,n)
  n=input%actual_matrix_dimension
  input_nnz=GetMatrixSize(input)
  IF(kind==0) THEN
    CALL CopyMatrix(input,a)
  ELSE
    CALL GetMatrixTripletList(input,real_entries);CALL CopyTripletList(real_entries,entries)
    IF(kind==2) THEN
      DO i=1,entries%CurrentSize
        theta=2.0_NTREAL*ACOS(-1.0_NTREAL)* &
          REAL(MOD(97*entries%DATA(i)%index_row,257)-MOD(97*entries%DATA(i)%index_column,257),NTREAL)/257
        entries%DATA(i)%point_value=entries%DATA(i)%point_value*CMPLX(COS(theta),SIN(theta),NTREAL)
      END DO
    END IF
    CALL ConstructEmptyMatrix(a,n,process_grid_in=grid,is_complex_in=.TRUE.)
    CALL FillMatrixFromTripletList(a,entries)
    CALL DestructTripletList(real_entries);CALL DestructTripletList(entries)
  END IF
  threshold=1e-8_NTREAL
  CALL MatrixMultiply(a,a,reference,threshold_in=threshold)
  ! Require the requested placement; a resident benchmark must never report CPU fallback as GPU time.
  CALL SelectGPUBackend(a,a,allowed,root_mode,1,input_nnz,preserve_padding=.TRUE.)
  IF(.NOT.allowed .OR. root_mode) CALL MPI_Abort(MPI_COMM_WORLD,10,ierr)
  CALL BeginGPUBackend(a,.FALSE.)
  CALL CheckGPUMemory(a,a,.FALSE.,allowed,input_nnz,preserve_padding=.TRUE.)
  IF(.NOT.allowed) CALL MPI_Abort(MPI_COMM_WORLD,11,ierr)
  ALLOCATE(seconds(reps),phases(4,reps),local_seconds(reps))
  IF(kind==0) THEN
    CALL UploadGPU(a,ga)
    ! Untimed readiness operation also establishes device residency for managed arrays.
    CALL MultiplyGPU(ga,ga,gresult,threshold)
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    DO rep=1,reps
      IF(profiling==1 .AND. rep==2) CALL profile_begin()
      t0=MPI_Wtime()
      CALL MultiplyGPU(ga,ga,gresult,threshold,phases(:,rep))
      seconds(rep)=MPI_Wtime()-t0
    END DO
    IF(profiling==1) CALL profile_end()
    local_seconds=phases(4,:)
    CALL ConstructEmptyMatrix(result,n,process_grid_in=grid)
    CALL DownloadGPU(gresult,result)
    CALL ReleaseGPU(ga);CALL ReleaseGPU(gresult)
  ELSE
    CALL GetMatrixTripletList(a,entries)
    nnz=entries%CurrentSize
    ALLOCATE(ir(MAX(1,nnz)),ic(MAX(1,nnz)),iv(MAX(1,nnz)))
    DO i=1,nnz
      ir(i)=entries%DATA(i)%index_row;ic(i)=entries%DATA(i)%index_column;iv(i)=entries%DATA(i)%point_value
    END DO
    CALL DestructTripletList(entries)
    rows=a%local_rows;columns=a%local_columns;dimension=a%logical_matrix_dimension
    row_offset=a%start_row-1;column_offset=a%start_column-1
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    status=resident_complex(rows,columns,dimension,row_offset,column_offset,nnz,ir,ic,iv, &
      threshold,reps,profiling,seconds,phases,output,count)
    IF(status/=0) CALL MPI_Abort(MPI_COMM_WORLD,12,ierr)
    local_seconds=phases(3,:)+phases(4,:)
    ALLOCATE(outrows(MAX(1,count)),outcols(MAX(1,count)),outvalues(MAX(1,count)))
    CALL complex_export(output,outrows,outcols,outvalues)
    CALL ConstructTripletList(entries,count)
    DO i=1,count
      entries%DATA(i)%index_row=outrows(i);entries%DATA(i)%index_column=outcols(i)
      entries%DATA(i)%point_value=outvalues(i)
    END DO
    CALL ConstructEmptyMatrix(result,n,process_grid_in=grid,is_complex_in=.TRUE.)
    CALL FillMatrixFromTripletList(result,entries,preduplicated_in=.TRUE.,prepartitioned_in=.TRUE.)
    CALL DestructTripletList(entries)
  END IF
  CALL EndGPUBackend()
  CALL MPI_Allreduce(MPI_IN_PLACE,seconds,reps,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,local_seconds,reps,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,phases,4*reps,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  output_nnz=GetMatrixSize(result)
  CALL CopyMatrix(result,difference);CALL IncrementMatrix(reference,difference,alpha_in=-1.0_NTREAL)
  error=MatrixNorm(difference)/MAX(1.0_NTREAL,MatrixNorm(reference))
  IF(.NOT.(error<1e-9_NTREAL)) CALL MPI_Abort(MPI_COMM_WORLD,13,ierr)
  IF(rank==0) THEN
    WRITE(*,'(A,5(1X,I0))') 'resident_input',n,nr,kind,input_nnz,output_nnz
    DO rep=1,reps
      WRITE(*,'(A,2(1X,I0),7(1X,ES24.16))') 'resident_product',kind,rep,seconds(rep),local_seconds(rep), &
        phases(:,rep),error
    END DO
    PRINT *, 'RESIDENT PRODUCT BENCH COMPLETE'
  END IF
  CALL DestructMatrix(input);CALL DestructMatrix(a);CALL DestructMatrix(reference)
  CALL DestructMatrix(result);CALL DestructMatrix(difference);CALL DestructProcessGrid(grid)
  CALL MPI_Finalize(ierr)
END PROGRAM
