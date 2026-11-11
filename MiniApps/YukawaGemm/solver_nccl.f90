MODULE SolverNCCL
  USE ISO_C_BINDING
  USE MPI
  IMPLICIT NONE
  LOGICAL :: nccl_active=.FALSE.
  INTERFACE
    INTEGER(C_INT) FUNCTION available() BIND(C,NAME='solver_nccl_available_')
      IMPORT C_INT
    END FUNCTION
    INTEGER(C_INT) FUNCTION get_id(id) BIND(C,NAME='solver_nccl_id_')
      IMPORT C_INT,C_INT8_T
      INTEGER(C_INT8_T) :: id(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION init(which,size,rank,id) BIND(C,NAME='solver_nccl_init_')
      IMPORT C_INT,C_INT8_T
      INTEGER(C_INT) :: which,size,rank
      INTEGER(C_INT8_T) :: id(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION set_members(which,members) BIND(C,NAME='solver_nccl_members_')
      IMPORT C_INT
      INTEGER(C_INT) :: which,members(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION gather_i(which,send,count,receive,counts,displs) BIND(C,NAME='solver_nccl_gather_i_')
      IMPORT C_INT
      INTEGER(C_INT) :: which,send(*),count,receive(*),counts(*),displs(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION gather_r(which,send,count,receive,counts,displs) BIND(C,NAME='solver_nccl_gather_r_')
      IMPORT C_INT,C_DOUBLE
      INTEGER(C_INT) :: which,count,counts(*),displs(*)
      REAL(C_DOUBLE) :: send(*),receive(*)
    END FUNCTION
    INTEGER(C_INT) FUNCTION sum_scalar(value) BIND(C,NAME='solver_nccl_sum_')
      IMPORT C_INT,C_DOUBLE
      REAL(C_DOUBLE) :: value
    END FUNCTION
    INTEGER(C_INT) FUNCTION finish() BIND(C,NAME='solver_nccl_finish_')
      IMPORT C_INT
    END FUNCTION
  END INTERFACE
CONTAINS
  SUBROUTINE CheckNCCL(status)
    INTEGER, INTENT(IN) :: status
    INTEGER :: ierr
    IF(status==0) RETURN
    PRINT *, 'NCCL communication failed; status:',status
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END SUBROUTINE

  SUBROUTINE InitSolverNCCL(row_comm,col_comm,world_comm)
    INTEGER, INTENT(IN) :: row_comm,col_comm
    INTEGER, OPTIONAL, INTENT(IN) :: world_comm
    INTEGER :: bootstrap_comm
    INTEGER :: shared,local_size,world_size,rank,ierr,which,comms(3),nr,cr
    INTEGER, ALLOCATABLE :: members(:)
    CHARACTER(16) :: shared_mode,multinode_mode
    INTEGER(C_INT8_T) :: id(128)
    IF(available()==0) CALL CheckNCCL(-1)
    bootstrap_comm=MPI_COMM_WORLD
    IF(PRESENT(world_comm)) bootstrap_comm=world_comm
    CALL MPI_Comm_size(bootstrap_comm,world_size,ierr)
    CALL MPI_Comm_rank(bootstrap_comm,rank,ierr)
    CALL MPI_Comm_split_type(bootstrap_comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,shared,ierr)
    CALL MPI_Comm_size(shared,local_size,ierr)
    CALL MPI_Comm_free(shared,ierr)
    CALL GET_ENVIRONMENT_VARIABLE('YUKAWA_NCCL_ALLOW_MULTINODE',multinode_mode)
    IF(world_size>4 .OR. local_size/=world_size) THEN
      IF(TRIM(multinode_mode)/='1') THEN
        IF(rank==0) PRINT *, 'Multi-node NCCL requires YUKAWA_NCCL_ALLOW_MULTINODE=1'
        CALL CheckNCCL(-2)
      END IF
      IF(rank==0) PRINT *, 'Experimental multi-node NCCL enabled; world/local ranks:',world_size,local_size
    END IF
    comms=[bootstrap_comm,row_comm,col_comm]
    CALL GET_ENVIRONMENT_VARIABLE('YUKAWA_NCCL_SHARED_COMM',shared_mode)
    DO which=0,2
      CALL MPI_Comm_size(comms(which+1),nr,ierr)
      CALL MPI_Comm_rank(comms(which+1),cr,ierr)
      id=0
      IF(cr==0 .AND. nr>1 .AND. (which==0 .OR. TRIM(shared_mode)/='1')) CALL CheckNCCL(get_id(id))
      CALL MPI_Bcast(id,SIZE(id),MPI_BYTE,0,comms(which+1),ierr)
      CALL CheckNCCL(init(which,nr,cr,id))
      ALLOCATE(members(nr))
      CALL MPI_Allgather(rank,1,MPI_INTEGER,members,1,MPI_INTEGER,comms(which+1),ierr)
      CALL CheckNCCL(set_members(which,members))
      DEALLOCATE(members)
    END DO
    nccl_active=.TRUE.
  END SUBROUTINE

  SUBROUTINE WarmSolverNCCL(row_comm,col_comm)
    INTEGER, INTENT(IN) :: row_comm,col_comm
    INTEGER :: which,nr,cr,i,ierr,comm
    INTEGER, ALLOCATABLE :: counts(:),displs(:),si(:),ri(:)
    REAL(C_DOUBLE), ALLOCATABLE :: sr(:),rr(:)
    REAL(C_DOUBLE) :: scalar
    scalar=0
    CALL SolverScalarSum(scalar)
    DO which=1,2
      comm=row_comm
      IF(which==2) comm=col_comm
      CALL MPI_Comm_size(comm,nr,ierr)
      CALL MPI_Comm_rank(comm,cr,ierr)
      ALLOCATE(counts(nr),displs(nr),si(nr),ri(nr*(nr+1)/2),sr(nr),rr(nr*(nr+1)/2))
      si=0; sr=0
      CALL PanelFixedInt(which,si,1,ri,comm)
      counts=1
      DO i=1,nr
        displs(i)=i-1
      END DO
      CALL PanelReal(which,sr,1,rr,counts,displs,comm)
      DO i=1,nr
        counts(i)=i
        displs(i)=i*(i-1)/2
      END DO
      CALL PanelInt(which,si,cr+1,ri,counts,displs,comm)
      CALL PanelReal(which,sr,cr+1,rr,counts,displs,comm)
      DEALLOCATE(counts,displs,si,ri,sr,rr)
    END DO
  END SUBROUTINE

  SUBROUTINE PanelFixedInt(which,send,count,receive,comm)
    INTEGER :: which,send(*),count,receive(*),comm
    INTEGER :: nr,i,ierr
    INTEGER, ALLOCATABLE :: counts(:),displs(:)
    IF(.NOT.nccl_active) THEN
      CALL MPI_Allgather(send,count,MPI_INTEGER,receive,count,MPI_INTEGER,comm,ierr)
      RETURN
    END IF
    CALL MPI_Comm_size(comm,nr,ierr)
    ALLOCATE(counts(nr),displs(nr))
    counts=count
    DO i=1,nr
      displs(i)=(i-1)*count
    END DO
    CALL CheckNCCL(gather_i(which,send,count,receive,counts,displs))
  END SUBROUTINE

  SUBROUTINE PanelInt(which,send,count,receive,counts,displs,comm)
    INTEGER :: which,send(*),count,receive(*),counts(*),displs(*),comm
    INTEGER :: ierr
    IF(nccl_active) THEN
      CALL CheckNCCL(gather_i(which,send,count,receive,counts,displs))
    ELSE
      CALL MPI_Allgatherv(send,count,MPI_INTEGER,receive,counts,displs,MPI_INTEGER,comm,ierr)
    END IF
  END SUBROUTINE

  SUBROUTINE PanelReal(which,send,count,receive,counts,displs,comm)
    INTEGER :: which,count,counts(*),displs(*),comm
    REAL(C_DOUBLE) :: send(*),receive(*)
    INTEGER :: ierr
    IF(nccl_active) THEN
      CALL CheckNCCL(gather_r(which,send,count,receive,counts,displs))
    ELSE
      CALL MPI_Allgatherv(send,count,MPI_DOUBLE_PRECISION,receive,counts,displs,MPI_DOUBLE_PRECISION,comm,ierr)
    END IF
  END SUBROUTINE

  SUBROUTINE SolverScalarSum(value)
    REAL(C_DOUBLE) :: value
    INTEGER :: ierr
    IF(nccl_active) THEN
      CALL CheckNCCL(sum_scalar(value))
    ELSE
      CALL MPI_Allreduce(MPI_IN_PLACE,value,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    END IF
  END SUBROUTINE

  SUBROUTINE FinishSolverNCCL()
    IF(nccl_active) CALL CheckNCCL(finish())
    nccl_active=.FALSE.
  END SUBROUTINE
END MODULE
