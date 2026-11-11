!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> YukawaGemm: a mini application for benchmarking NTPoly's distributed
!> sparse matrix-matrix multiply.
!>
!> It reads a synthetic sparse matrix (a Hamiltonian-like matrix generated
!> from a screened Coulomb (Yukawa) potential on a realistic atomic geometry
!> by generate_matrices.py in this directory), then repeatedly multiplies it
!> by itself under a sparsity threshold, timing every multiply. This mirrors
!> the inner loop of density matrix purification, whose cost is dominated by
!> a long chain of sparse multiplies, and is therefore the primary kernel to
!> track while porting NTPoly to new hardware.
!>
!> Two multiply modes are supported:
!>   square (default): B <- B*B each round, tracking the fill evolution the
!>                     way a purification iteration chain does.
!>   fixed:            C <- A*A each round, a steady state multiply with
!>                     constant sparsity.
!>
!> Every run writes a JSON summary (timings, nonzero counts, and a checksum
!> fingerprint) so results can be compared across processor counts and,
!> eventually, across CPU and GPU backends.
PROGRAM YukawaGemm
  USE DataTypesModule, ONLY : NTREAL, NTLONG
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixTrace, &
       & MatrixGrandSum, MeasureAsymmetry, ScaleMatrix, IncrementMatrix, &
       & DotMatrix, MatrixNorm, MatrixSigma
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE DensityMatrixSolversModule, ONLY : TRS2
  USE SolverParametersModule
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation, &
       & DestructPermutation
  USE LoadBalancerModule, ONLY : PermuteMatrix
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool
  !! Local (non-distributed) modules — used for the --direct path that
  !! bypasses the distributed wrapper at np=1, eliminating gather/compose
  !! overhead that thrashes managed-memory pages between CPU and GPU.
  USE SMatrixModule, ONLY: Matrix_lsr, &
       & LocalTranspose => TransposeMatrix, &
       & LocalConstructEmpty => ConstructEmptyMatrix, &
       & LocalCopy => CopyMatrix, &
       & LocalDestruct => DestructMatrix
  USE SMatrixAlgebraModule, ONLY: LocalGemm => MatrixMultiply
  USE MatrixMemoryPoolModule, ONLY: MatrixMemoryPool_lr, &
       & LocalConstructPool => ConstructMatrixMemoryPool, &
       & LocalDestructPool => DestructMatrixMemoryPool
  USE, INTRINSIC :: ISO_C_BINDING
  USE LoggingModule, ONLY: ActivateLogger, DeactivateLogger
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE SolverNCCL, ONLY: InitSolverNCCL, WarmSolverNCCL, FinishSolverNCCL, PanelFixedInt, PanelInt, PanelReal, SolverScalarSum
  USE SolverCSR
  USE ConvergenceMonitor
  USE openacc
  USE MPI
  IMPLICIT NONE

  !! cuSPARSE hybrid interface (C wrapper via ISO_C_BINDING).
  !! Computes C = A*A with cuSPARSE SpGEMM, then thresholds.
  INTERFACE
     SUBROUTINE gpu_profile_begin() BIND(C, NAME='gpu_profile_begin_')
     END SUBROUTINE
     SUBROUTINE gpu_profile_end() BIND(C, NAME='gpu_profile_end_')
     END SUBROUTINE
     INTEGER(C_INT) FUNCTION gpu_resident_square(input_slot,output_slot,m,n,k,pr,pc,thr,times) &
          BIND(C,NAME='gpu_resident_square_')
       IMPORT :: C_INT,C_DOUBLE
       INTEGER(C_INT), VALUE :: input_slot,output_slot,m,n,k,pr,pc
       REAL(C_DOUBLE), VALUE :: thr
       REAL(C_DOUBLE), INTENT(OUT) :: times(4)
     END FUNCTION
     INTEGER(C_INT) FUNCTION gpu_resident_free() BIND(C,NAME='gpu_resident_free_')
       IMPORT :: C_INT
     END FUNCTION
     INTEGER(C_INT) FUNCTION cusparse_solver_begin(m,k) BIND(C, NAME='cusparse_solver_begin_')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE :: m,k
     END FUNCTION
     INTEGER(C_INT) FUNCTION cusparse_gemm_threshold( &
          & n, nnz_a, row_ptr_a, col_ind_a, val_a, threshold, &
          & nnz_c, row_ptr_c, col_ind_c, val_c) &
          & BIND(C, NAME='cusparse_gemm_threshold_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT), VALUE :: n, nnz_a
       INTEGER(C_INT), INTENT(IN) :: row_ptr_a(*), col_ind_a(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_a(*)
       REAL(C_DOUBLE), VALUE :: threshold
       INTEGER(C_INT), INTENT(OUT) :: nnz_c
       TYPE(C_PTR), INTENT(OUT) :: row_ptr_c, col_ind_c, val_c
     END FUNCTION cusparse_gemm_threshold

     SUBROUTINE cusparse_gemm_free(row_ptr, col_ind, val) &
          & BIND(C, NAME='cusparse_gemm_free_')
       IMPORT :: C_PTR
       TYPE(C_PTR), VALUE :: row_ptr, col_ind, val
     END SUBROUTINE cusparse_gemm_free

     !! Multi-GPU path: C = A_local * B, then threshold.
     !! A_local is m×n, B is n×k, C is m×k.
     INTEGER(C_INT) FUNCTION cusparse_gemm_ab_threshold( &
          & m, n, k, nnz_a, row_ptr_a, col_ind_a, val_a, &
          & nnz_b, row_ptr_b, col_ind_b, val_b, &
          & threshold, &
          & nnz_c, row_ptr_c, col_ind_c, val_c) &
          & BIND(C, NAME='cusparse_gemm_ab_threshold_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT), VALUE :: m, n, k, nnz_a, nnz_b
       INTEGER(C_INT), INTENT(IN) :: row_ptr_a(*), col_ind_a(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_a(*)
       INTEGER(C_INT), INTENT(IN) :: row_ptr_b(*), col_ind_b(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_b(*)
       REAL(C_DOUBLE), VALUE :: threshold
       INTEGER(C_INT), INTENT(OUT) :: nnz_c
       TYPE(C_PTR), INTENT(OUT) :: row_ptr_c, col_ind_c, val_c
     END FUNCTION cusparse_gemm_ab_threshold

     INTEGER(C_INT) FUNCTION cusparse_gemm_ab_dense( &
          & m, n, k, nnz_a, row_ptr_a, col_ind_a, val_a, &
          & nnz_b, row_ptr_b, col_ind_b, val_b, &
          & threshold, &
          & nnz_c, row_ptr_c, col_ind_c, val_c) &
          & BIND(C, NAME='cusparse_gemm_ab_dense_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT), VALUE :: m, n, k, nnz_a, nnz_b
       INTEGER(C_INT), INTENT(IN) :: row_ptr_a(*), col_ind_a(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_a(*)
       INTEGER(C_INT), INTENT(IN) :: row_ptr_b(*), col_ind_b(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_b(*)
       REAL(C_DOUBLE), VALUE :: threshold
       INTEGER(C_INT), INTENT(OUT) :: nnz_c
       TYPE(C_PTR), INTENT(OUT) :: row_ptr_c, col_ind_c, val_c
     END FUNCTION cusparse_gemm_ab_dense

     !! 2.5D k-split path: partial product P = A_local * B (no threshold),
     !! downloaded to host persistent buffers (exchanged over MPI, then
     !! merged on GPU by cusparse_ksplit_merge).
     INTEGER(C_INT) FUNCTION cusparse_gemm_ab_partial( &
          & m, n, k, nnz_a, row_ptr_a, col_ind_a, val_a, &
          & nnz_b, row_ptr_b, col_ind_b, val_b, &
          & nnz_p, row_ptr_p, col_ind_p, val_p) &
          & BIND(C, NAME='cusparse_gemm_ab_partial_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT), VALUE :: m, n, k, nnz_a, nnz_b
       INTEGER(C_INT), INTENT(IN) :: row_ptr_a(*), col_ind_a(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_a(*)
       INTEGER(C_INT), INTENT(IN) :: row_ptr_b(*), col_ind_b(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_b(*)
       INTEGER(C_INT), INTENT(OUT) :: nnz_p
       TYPE(C_PTR), INTENT(OUT) :: row_ptr_p, col_ind_p, val_p
     END FUNCTION cusparse_gemm_ab_partial

     !! 2.5D k-split path: merge ps partial products (CSR, concatenated and
     !! offset-fixed on the host) on the GPU: dense scatter-add, threshold,
     !! compact.  Returns the thresholded merged block in persistent host
     !! buffers, same layout as cusparse_gemm_ab_threshold output.
     INTEGER(C_INT) FUNCTION cusparse_ksplit_merge( &
          & ps, m, kk, total_nnz, row_ptr_all, col_ind_all, val_all, &
          & threshold, nnz_c, row_ptr_c, col_ind_c, val_c) &
          & BIND(C, NAME='cusparse_ksplit_merge_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT), VALUE :: ps, m, kk, total_nnz
       INTEGER(C_INT), INTENT(IN) :: row_ptr_all(*), col_ind_all(*)
       REAL(C_DOUBLE), INTENT(IN) :: val_all(*)
       REAL(C_DOUBLE), VALUE :: threshold
       INTEGER(C_INT), INTENT(OUT) :: nnz_c
       TYPE(C_PTR), INTENT(OUT) :: row_ptr_c, col_ind_c, val_c
     END FUNCTION cusparse_ksplit_merge

     !! P1 GPU-resident combine ops (managed-memory slots, CSR, rows sorted).
     SUBROUTINE gpu_ops_init(err, dev_out) &
          & BIND(C, NAME='gpu_ops_init_')
       IMPORT :: C_INT
       INTEGER(C_INT) :: err, dev_out
     END SUBROUTINE gpu_ops_init

     SUBROUTINE gpu_ops_alloc(slot, m, nnz, rp_out, ci_out, val_out, err) &
          & BIND(C, NAME='gpu_ops_alloc_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT) :: slot, m, nnz, err
       TYPE(C_PTR) :: rp_out, ci_out, val_out
     END SUBROUTINE gpu_ops_alloc

     SUBROUTINE gpu_ops_scale(slot, alpha) &
          & BIND(C, NAME='gpu_ops_scale_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot
       REAL(C_DOUBLE) :: alpha
     END SUBROUTINE gpu_ops_scale

     SUBROUTINE gpu_ops_increment(slot_c, slot_b, alpha, nnz_out) &
          & BIND(C, NAME='gpu_ops_increment_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot_c, slot_b, nnz_out
       REAL(C_DOUBLE) :: alpha
     END SUBROUTINE gpu_ops_increment

     SUBROUTINE gpu_ops_increment_threshold(slot_c, slot_b, alpha, thr, nnz_out) &
          & BIND(C, NAME='gpu_ops_increment_threshold_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot_c, slot_b, nnz_out
       REAL(C_DOUBLE) :: alpha, thr
     END SUBROUTINE

     SUBROUTINE gpu_ops_copy(dst, src) &
          & BIND(C, NAME='gpu_ops_copy_')
       IMPORT :: C_INT
       INTEGER(C_INT) :: dst, src
     END SUBROUTINE gpu_ops_copy

     SUBROUTINE gpu_ops_grandsum(slot, sum_out) &
          & BIND(C, NAME='gpu_ops_grandsum_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot
       REAL(C_DOUBLE) :: sum_out
     END SUBROUTINE gpu_ops_grandsum

     SUBROUTINE gpu_ops_trace(slot, diag_offset, trace_out) &
          & BIND(C, NAME='gpu_ops_trace_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot, diag_offset
       REAL(C_DOUBLE) :: trace_out
     END SUBROUTINE gpu_ops_trace

     SUBROUTINE gpu_ops_dot(slot_a, slot_b, dot_out) &
          & BIND(C, NAME='gpu_ops_dot_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot_a, slot_b
       REAL(C_DOUBLE) :: dot_out
     END SUBROUTINE gpu_ops_dot

     SUBROUTINE gpu_ops_colsum(slot, kk, colsum_out) &
          & BIND(C, NAME='gpu_ops_colsum_')
       IMPORT :: C_INT, C_DOUBLE
       INTEGER(C_INT) :: slot, kk
       REAL(C_DOUBLE) :: colsum_out(*)
     END SUBROUTINE gpu_ops_colsum

     SUBROUTINE gpu_ops_get(slot, rp, ci, val, m, nnz) &
          & BIND(C, NAME='gpu_ops_get_')
       IMPORT :: C_INT, C_DOUBLE, C_PTR
       INTEGER(C_INT) :: slot
       TYPE(C_PTR) :: rp, ci, val
       INTEGER(C_INT) :: m, nnz
     END SUBROUTINE gpu_ops_get
  END INTERFACE
  !! Parameters
  CHARACTER(LEN = 1024) :: input_file, output_file, mode_string
  INTEGER :: num_rounds, num_warmup
  REAL(NTREAL) :: threshold
  INTEGER :: process_rows, process_columns, process_slices
  LOGICAL :: have_grid, do_convert
  !! Matrices
  TYPE(Matrix_ps) :: matA, matB, matC
  TYPE(MatrixMemoryPool_p) :: memory_pool
  !! Direct-path local matrices (bypass distributed wrapper at np=1)
  TYPE(Matrix_lsr) :: matAT_local, matC_local
  TYPE(MatrixMemoryPool_lr) :: local_pool
  LOGICAL :: use_direct
  LOGICAL :: use_cusparse
  LOGICAL :: use_cusparse_mpi
  REAL(NTREAL) :: sparsity_est
  !! cuSPARSE hybrid path variables
  INTEGER(C_INT) :: cs_status, cs_nnz_c, cs_nnz_a
  TYPE(C_PTR) :: cs_row_ptr_c, cs_col_ind_c, cs_val_c
  INTEGER(C_INT), POINTER :: cs_f_row_ptr(:), cs_f_col_ind(:)
  REAL(C_DOUBLE), POINTER :: cs_f_val(:)
  LOGICAL :: cs_ok
  !! cuSPARSE multi-GPU path variables (1D row decomposition)
  INTEGER :: n_full, m_loc, start_row, end_row
  INTEGER(NTLONG) :: nnz_full, nnz_loc
  INTEGER, ALLOCATABLE :: row_ptr_full(:), col_ind_full(:)
  REAL(NTREAL), ALLOCATABLE :: val_full(:)
  INTEGER, ALLOCATABLE :: row_ptr_loc(:), col_ind_loc(:)
  REAL(NTREAL), ALLOCATABLE :: val_loc(:)
  INTEGER(C_INT) :: cs_ab_status, cs_ab_nnz_c, cs_ab_nnz_a, cs_ab_nnz_b
  TYPE(C_PTR) :: cs_ab_row_ptr_c, cs_ab_col_ind_c, cs_ab_val_c
  INTEGER(C_INT), POINTER :: cs_ab_f_row_ptr(:), cs_ab_f_col_ind(:)
  REAL(C_DOUBLE), POINTER :: cs_ab_f_val(:)
  LOGICAL :: use_symmetric
  LOGICAL :: use_load_balance
  INTEGER :: partner_rank
  TYPE(Permutation_t) :: balance_perm
  TYPE(Matrix_ps) :: matA_orig
  INTEGER, ALLOCATABLE :: seed_array(:)
  INTEGER :: seed_size
  INTEGER(C_INT), ALLOCATABLE, TARGET :: ct_row_ptr(:), ct_col_ind(:)
  REAL(C_DOUBLE), ALLOCATABLE, TARGET :: ct_vals(:)
  INTEGER(C_INT), ALLOCATABLE, TARGET :: recv_row_ptr(:), recv_col_ind(:)
  REAL(C_DOUBLE), ALLOCATABLE, TARGET :: recv_vals(:)
  INTEGER :: ct_nnz, recv_nnz
  INTEGER, ALLOCATABLE :: ct_counts(:), ct_offset(:)
  REAL(8) :: tpre, tspgemm, ttranspose, tsend, trecv
  REAL(NTREAL) :: local_trace, local_grand_sum, global_trace, global_grand_sum
  INTEGER(NTLONG) :: local_nnz, global_nnz
  INTEGER :: lr_loc, kk
  !! cuSPARSE 2D SUMMA path variables
  INTEGER :: pr_loc, pc_loc, ps_loc      !! process grid dims
  INTEGER :: my_row, my_col             !! position in grid
  INTEGER :: row_comm, col_comm         !! MPI communicators
  INTEGER :: within_comm                !! within-slice comm (fingerprint reduce)
  INTEGER :: JJ                          !! loop variable
  INTEGER :: offset_acc, idx_base       !! outer_index offset fixup
  INTEGER :: kk2, global_row             !! fingerprint loop vars
  INTEGER :: n_row_blk, n_col_blk        !! row-block dims (N/pr × N)
  INTEGER :: row_nnz_loc, row_nnz_total  !! row-block nnz counts
  INTEGER :: col_row_nnz_total           !! nnz of row-block from row my_col
  INTEGER :: my_col_offset               !! global column offset of this column block
  INTEGER :: r_pos, c_pos                !! transpose loop vars
  INTEGER :: mpi_status_buf(MPI_STATUS_SIZE) !! MPI status for Sendrecv
  !! Gathered row-block (my row) in CSC: m_loc × N
  INTEGER, ALLOCATABLE, TARGET :: row_blk_outer(:), row_blk_inner(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: row_blk_vals(:)
  !! My full local block (all sub-blocks merged) in CSC: m_loc × kk
  INTEGER, ALLOCATABLE, TARGET :: loc_outer(:), loc_inner(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: loc_vals(:)
  !! Row-block from row my_col (received via Sendrecv) in CSC: m_loc × N
  !! For symmetric A, this CSC = B's CSR (N × m_loc), no transpose needed.
  INTEGER, ALLOCATABLE, TARGET :: col_row_blk_outer(:), col_row_blk_inner(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: col_row_blk_vals(:)
  !! A_local = transpose of row-block from CSR→CSC, in CSR: m_loc × N
  INTEGER, ALLOCATABLE :: aloc_row_ptr(:), aloc_col_ind(:)
  REAL(NTREAL), ALLOCATABLE :: aloc_vals(:)
  INTEGER, ALLOCATABLE :: row_count(:)   !! scratch for transpose
  !! B = row-block from my_col in CSC, aliased as CSR (N × m_loc)
  INTEGER, POINTER :: b_row_ptr(:), b_col_ind(:)
  REAL(NTREAL), POINTER :: b_vals(:)
  !! Row-block AllGather helpers (along row_comm)
  INTEGER, ALLOCATABLE :: row_nnz_counts(:), row_displs(:)
  !! 2.5D k-split variables (--ksplit): each slice s computes the partial
  !! C_ij^(s) over the s-th slab of the shared dimension, partials are
  !! exchanged over between_slice_comm and merged on the GPU.
  LOGICAL :: use_ksplit
  INTEGER :: my_slice, btw_comm, btw_rank, btw_size
  INTEGER :: kfull, ks_lo, ks_b0, ks_b1, ks_pos, ks_base
  INTEGER :: nnz_a_s, nnz_b_s, ks_nnz_p, ks_total
  INTEGER(C_INT) :: ks_nnz_p_c
  TYPE(C_PTR) :: ks_rp_p, ks_ci_p, ks_va_p
  INTEGER(C_INT), POINTER :: ks_f_rp(:), ks_f_ci(:)
  REAL(C_DOUBLE), POINTER :: ks_f_val(:)
  INTEGER, ALLOCATABLE, TARGET :: aloc_s_row_ptr(:), aloc_s_col_ind(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: aloc_s_vals(:)
  INTEGER, ALLOCATABLE, TARGET :: b_s_row_ptr(:), b_s_col_ind(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: b_s_vals(:)
  INTEGER, ALLOCATABLE, TARGET :: ksp_row_ptr(:), ksp_col_ind(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: ksp_vals(:)
  INTEGER, ALLOCATABLE, TARGET :: ks_my_rp(:)
  REAL(NTREAL), ALLOCATABLE, TARGET :: ks_dummy_val(:)
  INTEGER, ALLOCATABLE :: ks_counts(:), ks_displs(:)
  INTEGER, ALLOCATABLE :: ks_rp_counts(:), ks_rp_displs(:)
  REAL(8) :: tpart, texch, tmerge
  REAL(8) :: tstage_read, tstage_gather, tstage_transpose, tstage_bcast
  CHARACTER(LEN = 32) :: backend_str
  !! P1: GPU-resident combine-op benchmarks (--op)
  LOGICAL :: use_op_bench
  CHARACTER(LEN = 32) :: op_mode
  INTEGER :: op_reps, op_warmup
  REAL(8) :: op_alpha
  INTEGER(C_INT) :: op_nnz_out
  INTEGER(C_INT) :: ops_slot1, ops_slot2, ops_err
  TYPE(C_PTR) :: op_rp_ptr, op_ci_ptr, op_val_ptr
  INTEGER(C_INT), POINTER :: mgmt_row_ptr(:), mgmt_col_ind(:)
  REAL(C_DOUBLE), POINTER :: mgmt_vals(:)
  REAL(C_DOUBLE) :: op_gs0, op_gs1
  REAL(C_DOUBLE) :: op_scalar_gpu, op_scalar_cpu
  REAL(NTREAL), ALLOCATABLE :: colsum(:)
  INTEGER :: diag_offset
  !! P3: TRS2 solver benchmark (--solver trs2)
  LOGICAL :: use_nccl, use_device_partials, warm_nccl, use_gpu_dense, use_resident_panels
  REAL(C_DOUBLE) :: resident_phase(4)
  REAL(8) :: nccl_setup_seconds, solver_setup_seconds, solver_setup_t0
  REAL(8) :: gpu_bind_seconds, gpu_restore_seconds, gpu_chain_seconds
  INTEGER :: solver_repetitions=1, solver_rep
  INTEGER(C_INT), PARAMETER :: solver_slot_initial=4
  CHARACTER(LEN=1024) :: solver_profile_path
  CHARACTER(LEN=32) :: solver_rep_string
  LOGICAL :: use_solver
  INTEGER :: max_iter, siter
  REAL(NTREAL) :: trace_target, conv_diff
  REAL(NTREAL) :: e_min, e_max, tr_val, energy_new, energy_old
  REAL(NTREAL) :: sigma_sign
  TYPE(Matrix_ps) :: matWH, matX, matX2, ISQ_mat, K_out
  TYPE(Matrix_lsr) :: solver_local_x
  TYPE(MatrixMemoryPool_p) :: solver_pool
  TYPE(SolverParameters_t) :: solver_params
  REAL(NTREAL) :: solver_energy
  INTEGER(C_INT) :: solver_slot_wh, solver_slot_x, solver_slot_x2
  REAL(C_DOUBLE) :: alpha_tmp
  INTEGER(C_INT) :: op_m_out
  INTEGER :: lr_local, ops_dev
  INTEGER :: total_iters
  CHARACTER(LEN = 32) :: solver_mode
  REAL(8) :: t_mult, t_red, t_upd, t_io
  REAL(8) :: phase(17), phase_max(17), profile_t0, iteration_t0, loop_t0
  REAL(8) :: gpu_solver_seconds, ref_solver_seconds
  REAL(8), ALLOCATABLE :: phase_history(:,:), all_phase_history(:,:,:)
  INTEGER :: profile_rank, profile_iteration
  INTEGER :: profile_unit
  REAL(8), ALLOCATABLE :: iter_times(:)
  !! P3 chunked multiply state
  TYPE(Monitor_t) :: gpu_monitor
  LOGICAL :: chunk_ok, gpu_converged
  REAL(NTREAL) :: reference_trace, reference_energy, parity_energy, parity_trace
  INTEGER :: log_ios
  CHARACTER(LEN=1024) :: solver_log_path
  INTEGER(C_INT64_T) :: chunk_work_budget = 4000000_C_INT64_T
  REAL(NTREAL) :: reference_norm, parity_matrix
  INTEGER, ALLOCATABLE :: chunk_ends(:), chunk_offsets(:)
  INTEGER :: nchunk, kfull_c, ks_hi, base
  INTEGER :: nnz_ac, nnz_bc
  INTEGER, ALLOCATABLE :: ac_rp(:), ac_ci(:)
  REAL(NTREAL), ALLOCATABLE :: ac_val(:)
  INTEGER, ALLOCATABLE :: ck_rp(:), ck_ci(:)
  REAL(NTREAL), ALLOCATABLE :: ck_vals(:)
  INTEGER, ALLOCATABLE :: bc_rp_c(:), bc_ci_c(:), ac_ci_c(:)
  REAL(NTREAL), ALLOCATABLE :: bc_val_c(:), ac_val_c(:)
  REAL(8), ALLOCATABLE :: op_times(:)
  REAL(8) :: dt_op, gpu_med, gpu_mean, cpu_med, rel_err
  REAL(NTREAL) :: gs_gpu_global, gs_cpu_global
  INTEGER(NTLONG) :: op_nnz_start, op_nnz_end
  TYPE(Matrix_ps) :: matB_ps
  !! Bookkeeping
  INTEGER :: my_rank, comm_size, ierr, provided
  INTEGER :: II, NR, unit_handle, consumed
  LOGICAL :: is_root, is_square, timed
  REAL(8) :: t0, t1, dt_local, dt_gemm, dt_rot
  REAL(8), ALLOCATABLE :: gemm_times(:), rot_times(:)
  INTEGER(NTLONG), ALLOCATABLE :: round_nnz(:)
  INTEGER(NTLONG) :: nnz_input, nnz_current
  REAL(NTREAL) :: trace_value, grand_sum, asymmetry
  REAL(NTREAL) :: mean_time, median_time, min_time, max_time
  INTEGER :: matrix_dim
  CHARACTER(LEN = 1024) :: argument, argument_value, binary_target
  CHARACTER(LEN = 8) :: date_str
  CHARACTER(LEN = 10) :: time_str
  CHARACTER(LEN = 5) :: zone_str
  CHARACTER(LEN = 4096) :: line, numbuf

  !! Defaults
  input_file = ""
  output_file = ""
  mode_string = "fixed"
  num_rounds = 10
  num_warmup = 3
  threshold = 1.0_NTREAL * 1e-4
  have_grid = .FALSE.
  do_convert = .FALSE.
  use_direct = .FALSE.
  use_cusparse = .FALSE.
  use_cusparse_mpi = .FALSE.
  use_symmetric = .FALSE.
  use_load_balance = .FALSE.
  use_ksplit = .FALSE.
  use_op_bench = .FALSE.
  use_nccl = .FALSE.
  use_device_partials = .FALSE.
  warm_nccl = .FALSE.
  use_gpu_dense = .FALSE.
  use_resident_panels = .FALSE.
  nccl_setup_seconds = 0
  use_solver = .FALSE.
  solver_mode = "trs2"
  max_iter = 60
  conv_diff = 1.0e-6_NTREAL
  trace_target = 0.0_NTREAL
  solver_slot_wh = 1
  solver_slot_x = 2
  solver_slot_x2 = 3
  op_mode = "scale"
  op_reps = 20
  op_warmup = 3
  op_alpha = 1.0_8 + 1.0e-2_8
  op_nnz_out = 0
  ops_slot1 = 1
  ops_slot2 = 2
  partner_rank = -1
  cs_ok = .FALSE.

  !! Parse the command line arguments.
  II = 1
  DO WHILE (II .LE. COMMAND_ARGUMENT_COUNT())
     CALL GET_COMMAND_ARGUMENT(II, argument)
     consumed = 1
     SELECT CASE (TRIM(argument))
     CASE ('--input')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        input_file = TRIM(argument_value)
        consumed = 2
     CASE ('--out')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        output_file = TRIM(argument_value)
        consumed = 2
     CASE ('--rounds')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) num_rounds
        consumed = 2
     CASE ('--warmup')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) num_warmup
        consumed = 2
     CASE ('--threshold')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) threshold
        consumed = 2
     CASE ('--mode')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        mode_string = TRIM(argument_value)
        consumed = 2
     CASE ('--prows')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) process_rows
        have_grid = .TRUE.
        consumed = 2
     CASE ('--pcols')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) process_columns
        have_grid = .TRUE.
        consumed = 2
     CASE ('--pslices')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) process_slices
        have_grid = .TRUE.
        consumed = 2
     CASE ('--convert')
        do_convert = .TRUE.
     CASE ('--direct')
        use_direct = .TRUE.
     CASE ('--cusparse')
        use_cusparse = .TRUE.
     CASE ('--cusparse-mpi')
        use_cusparse_mpi = .TRUE.
     CASE ('--cusparse-2d')
        use_cusparse_mpi = .TRUE.
     CASE ('--symmetric')
        use_symmetric = .TRUE.
     CASE ('--ksplit')
        use_ksplit = .TRUE.
     CASE ('--op')
        use_op_bench = .TRUE.
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        op_mode = TRIM(argument_value)
        consumed = 2
     CASE ('--reps')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) op_reps
        consumed = 2
     CASE ('--alpha')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) op_alpha
        consumed = 2
     CASE ('--nccl')
        use_nccl = .TRUE.
     CASE ('--device-partials')
        use_device_partials = .TRUE.
     CASE ('--solver-reps')
        CALL GET_COMMAND_ARGUMENT(II+1,argument_value)
        READ(argument_value,*) solver_repetitions
        consumed=2
     CASE ('--resident-panels')
        use_resident_panels = .TRUE.
     CASE ('--gpu-dense')
        use_gpu_dense = .TRUE.
     CASE ('--nccl-warmup')
        warm_nccl = .TRUE.
     CASE ('--solver')
        use_solver = .TRUE.
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        solver_mode = TRIM(argument_value)
        consumed = 2
     CASE ('--chunk-work')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) chunk_work_budget
        consumed = 2
     CASE ('--nelec')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) trace_target
        consumed = 2
     CASE ('--max-iter')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) max_iter
        consumed = 2
     CASE ('--conv-diff')
        CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
        READ (argument_value, *) conv_diff
        consumed = 2
     CASE ('--load-balance')
        use_load_balance = .TRUE.
     CASE DEFAULT
        ! Ignore unknown arguments.
     END SELECT
     II = II + consumed
  END DO

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)
  is_root = my_rank .EQ. 0

  IF (TRIM(input_file) .EQ. "") THEN
     IF (is_root) THEN
        CALL PrintUsage()
     END IF
     CALL MPI_Finalize(ierr)
     STOP
  END IF

  IF (use_ksplit .AND. .NOT. use_cusparse_mpi) THEN
     IF (is_root) WRITE (*,*) "Error: --ksplit requires --cusparse-2d"
     CALL MPI_Finalize(ierr)
     STOP
  END IF

  !! Setup the process grid.
  IF (have_grid) THEN
     CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, &
          & process_columns, process_slices)
  ELSE
     CALL ConstructProcessGrid(MPI_COMM_WORLD)
  END IF

  !! Optional conversion mode: matrix market -> NTPoly binary.
  !! The binary format is much faster to read for large matrices.
  IF (do_convert) THEN
     CALL ReadInputMatrix(matA, TRIM(input_file))
     binary_target = TRIM(output_file)
     IF (TRIM(binary_target) .EQ. "") THEN
        binary_target = TRIM(input_file)
        IF (INDEX(TRIM(binary_target), ".", BACK = .TRUE.) .GT. 0) THEN
           binary_target = binary_target(1:INDEX(TRIM(binary_target), ".", &
                & BACK = .TRUE.) - 1)
        END IF
        binary_target = TRIM(binary_target) // ".bin"
     END IF
     CALL WriteMatrixToBinary(matA, TRIM(binary_target))
     IF (is_root) THEN
        WRITE (*, '(A,1X,A)') "Converted", TRIM(input_file)
        WRITE (*, '(A,1X,A)') "to     :", TRIM(binary_target)
        WRITE (*, '(A,1X,I0)') "nnz    :", GetMatrixSize(matA)
     END IF
     CALL DestructMatrix(matA)
     CALL DestructProcessGrid
     CALL MPI_Finalize(ierr)
     STOP
  END IF

   !! cuSPARSE 2D SUMMA path (requires square process grid pr = pc).
   !! For C = A*A with symmetric A:
   !!   1. AllGather row-block along row_comm -> A_rows (m_loc x N) in CSC
   !!   2. Sendrecv row-block along col_comm to get B's source (row my_col)
   !!   3. Transpose A_rows from CSC -> CSR -> A_local (m_loc x N)
   !!   4. B = row-block from my_col in CSC, passed as CSR (N x m_loc), no
   !!      transpose (CSC of m_loc x N = CSR of N x m_loc for the transpose,
   !!      since A is symmetric)
   !!   5. cuSPARSE: C = A_local * B (m_loc x m_loc), then threshold
   !! Output stays in 2D block layout -- no redistribution for purification.
   is_square = TRIM(mode_string) .NE. "fixed"
   IF (use_nccl .AND. .NOT. use_solver) CALL SolverAbort('--nccl requires --solver trs2')
   IF (use_device_partials .AND. .NOT. use_solver) CALL SolverAbort('--device-partials requires --solver trs2')
   IF (solver_repetitions<1 .OR. (solver_repetitions/=1 .AND. .NOT. use_solver)) &
        CALL SolverAbort('--solver-reps requires a solver and a positive count')
   IF (use_resident_panels .AND. (.NOT. use_gpu_dense .OR. .NOT. use_nccl)) &
        CALL SolverAbort('--resident-panels requires --gpu-dense and --nccl')
   IF (use_gpu_dense .AND. .NOT. use_solver) CALL SolverAbort('--gpu-dense requires --solver trs2')
   IF (warm_nccl .AND. .NOT. use_nccl) CALL SolverAbort('--nccl-warmup requires --nccl')
   IF (use_solver .AND. use_op_bench) CALL SolverAbort('choose --solver or --op')

   !! ---------------------------------------------------------------------
   !! P1 (GPU_ROADMAP.md): GPU-resident combine-op benchmark (--op).
   !! Reads the (optionally load-balanced) matrix in the natural 2D block
   !! layout, uploads each rank's local block to a device slot, then times
   !! N repetitions of scale / increment / copy entirely on device. The
   !! CPU reference times the same repetitions through NTPoly's
   !! ScaleMatrix / IncrementMatrix / CopyMatrix. Grand-sum fingerprints
   !! verify the entry set (scale: x alpha^reps, increment: x (1+alpha)^reps)
   !! and CPU-vs-GPU parity of the global grand sum.
   IF (use_op_bench .AND. .NOT. is_square) THEN
      IF (TRIM(op_mode) .NE. "scale" .AND. TRIM(op_mode) .NE. "increment" &
           & .AND. TRIM(op_mode) .NE. "copy" .AND. TRIM(op_mode) .NE. "trace" &
           & .AND. TRIM(op_mode) .NE. "dot" .AND. TRIM(op_mode) .NE. "norm" &
           & .AND. TRIM(op_mode) .NE. "sigma" .AND. TRIM(op_mode) .NE. "grandsum") THEN
         IF (is_root) WRITE (*,*) &
              & "Error: --op must be scale|increment|copy|trace|dot|norm|sigma|grandsum"
         CALL DestructProcessGrid
         CALL MPI_Finalize(ierr)
         STOP
      END IF

      !! Read (optionally load-balanced) matrix.
      t0 = MPI_Wtime()
      IF (use_load_balance) THEN
         CALL ReadInputMatrix(matA_orig, TRIM(input_file))
         matrix_dim = matA_orig%actual_matrix_dimension
         CALL RANDOM_SEED(SIZE = seed_size)
         ALLOCATE (seed_array(seed_size))
         DO II = 1, seed_size
            seed_array(II) = 424242 + II
         END DO
         CALL RANDOM_SEED(PUT = seed_array)
         DEALLOCATE (seed_array)
         CALL ConstructRandomPermutation(balance_perm, matrix_dim)
         CALL PermuteMatrix(matA_orig, matA, balance_perm)
         CALL DestructPermutation(balance_perm)
         CALL DestructMatrix(matA_orig)
      ELSE
         CALL ReadInputMatrix(matA, TRIM(input_file))
      END IF
      tstage_read = MPI_Wtime() - t0
      matrix_dim = matA%actual_matrix_dimension
      nnz_input = GetMatrixSize(matA)
      m_loc = matA%local_rows
      kk = matA%local_columns

      !! Merge sub-blocks into one CSC block (m_loc x kk), then transpose
      !! to CSR with 1-based LOCAL columns. The CSR arrays live in CUDA
      !! MANAGED memory: alloc first, C_F_POINTER over them, then build
      !! in place (CPU writes; GPU kernels fault the pages over later).
      !! Row-major emission order makes every row column-sorted by
      !! construction (the merge invariant).
      !! NOTE: gpu_ops_init (cudaSetDevice) must run AFTER PermuteMatrix —
      !! NTPoly's multiply is OpenACC-enabled and its runtime binds lazily;
      !! an early cudaSetDevice races it (CUDA_ERROR_INVALID_HANDLE).
      CALL gpu_ops_init(ops_err, ops_dev)
      IF (ops_err .NE. 0) THEN
         WRITE (*,*) "Error: gpu_ops_init (device bind) failed on rank", &
              & my_rank
         CALL DestructMatrix(matA)
         CALL DestructProcessGrid
         CALL MPI_Finalize(ierr)
         STOP
      END IF
      CALL MergeLocalBlocksToCSC(matA, loc_outer, loc_inner, loc_vals, &
           & row_nnz_loc)
      CALL gpu_ops_alloc(ops_slot1, m_loc, row_nnz_loc, op_rp_ptr, &
           & op_ci_ptr, op_val_ptr, ops_err)
      IF (ops_err .NE. 0) THEN
         WRITE (*,*) "Error: gpu_ops_alloc failed on rank", my_rank
         CALL DestructMatrix(matA)
         CALL DestructProcessGrid
         CALL MPI_Finalize(ierr)
         STOP
      END IF
      CALL C_F_POINTER(op_rp_ptr, mgmt_row_ptr, [m_loc + 1])
      CALL C_F_POINTER(op_ci_ptr, mgmt_col_ind, [MAX(row_nnz_loc, 1)])
      CALL C_F_POINTER(op_val_ptr, mgmt_vals, [MAX(row_nnz_loc, 1)])
      mgmt_row_ptr = 0
      DO II = 1, row_nnz_loc
         mgmt_row_ptr(loc_inner(II) + 1) = &
              & mgmt_row_ptr(loc_inner(II) + 1) + 1
      END DO
      mgmt_row_ptr(1) = 0
      DO II = 1, m_loc
         mgmt_row_ptr(II + 1) = mgmt_row_ptr(II) + mgmt_row_ptr(II + 1)
      END DO
      ALLOCATE (row_count(m_loc))
      row_count = 0
      DO JJ = 1, kk
         DO II = loc_outer(JJ) + 1, loc_outer(JJ + 1)
            r_pos = loc_inner(II)
            c_pos = mgmt_row_ptr(r_pos) + row_count(r_pos) + 1
            mgmt_col_ind(c_pos) = INT(JJ, KIND = C_INT)
            mgmt_vals(c_pos) = REAL(loc_vals(II), KIND = C_DOUBLE)
            row_count(r_pos) = row_count(r_pos) + 1
         END DO
      END DO
      DEALLOCATE (row_count)
      DEALLOCATE (loc_outer, loc_inner, loc_vals)

      !! For increment/copy/dot: slot 2 = device copy of X, plus an
      !! NTPoly copy (matB_ps) for the CPU reference path.
      IF (TRIM(op_mode) .EQ. "increment" .OR. TRIM(op_mode) .EQ. "copy" &
           & .OR. TRIM(op_mode) .EQ. "dot") THEN
         CALL gpu_ops_copy(ops_slot2, ops_slot1)
         CALL CopyMatrix(matA, matB_ps)
      END IF
      !! Reduction-mode setup: grid communicators for the MPI combines and
      !! the block-local diagonal offset for trace.
      my_row = matA%process_grid%my_row
      my_col = matA%process_grid%my_column
      col_comm = matA%process_grid%column_comm
      row_comm = matA%process_grid%row_comm
      diag_offset = (my_row - my_col) * m_loc
      ALLOCATE (colsum(MAX(kk, 1)))
      CALL gpu_ops_grandsum(ops_slot1, op_gs0)
      gs_gpu_global = 0.0_NTREAL
      CALL MPI_Allreduce(REAL(op_gs0, KIND = NTREAL), gs_gpu_global, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      op_nnz_start = row_nnz_loc

      !! Print parameters.
      IF (is_root) THEN
         WRITE (*, '(A)') "=== YukawaGemm (op bench) ==="
         WRITE (*, '(A,1X,A)') "input     :", TRIM(input_file)
         WRITE (*, '(A,1X,I0)') "dimension :", matrix_dim
         WRITE (*, '(A,1X,I0)') "input nnz :", nnz_input
         WRITE (*, '(A,1X,I0)') "mpi ranks :", comm_size
         WRITE (*, '(A,1X,A)') "op        :", TRIM(op_mode)
         WRITE (*, '(A,1X,I0)') "reps      :", op_reps
         WRITE (*, '(A,1X,ES10.3)') "alpha     :", op_alpha
         WRITE (*, '(A,1X,A)') "backend   : GPU resident combine ops"
         WRITE (*, '(A)') "-----------------"
      END IF

      !! GPU benchmark: warmup + timed reps, max-over-ranks per rep.
      ALLOCATE (op_times(op_reps))
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      DO II = 1, op_warmup + op_reps
         timed = II .GT. op_warmup
         NR = II - op_warmup
         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         t0 = MPI_Wtime()
         SELECT CASE (TRIM(op_mode))
         CASE ("scale")
            CALL gpu_ops_scale(ops_slot1, op_alpha)
         CASE ("increment")
            CALL gpu_ops_increment(ops_slot1, ops_slot2, op_alpha, op_nnz_out)
         CASE ("copy")
            CALL gpu_ops_copy(ops_slot2, ops_slot1)
         CASE ("trace")
            CALL gpu_ops_trace(ops_slot1, diag_offset, op_scalar_gpu)
            CALL MPI_Allreduce(MPI_IN_PLACE, op_scalar_gpu, 1, &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         CASE ("dot")
            CALL gpu_ops_dot(ops_slot1, ops_slot2, op_scalar_gpu)
            CALL MPI_Allreduce(MPI_IN_PLACE, op_scalar_gpu, 1, &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         CASE ("norm", "sigma")
            CALL gpu_ops_colsum(ops_slot1, kk, colsum)
            CALL MPI_Allreduce(MPI_IN_PLACE, colsum, kk, &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, col_comm, ierr)
            op_scalar_gpu = MAXVAL(colsum(1:kk))
            CALL MPI_Allreduce(MPI_IN_PLACE, op_scalar_gpu, 1, &
                 & MPI_DOUBLE_PRECISION, MPI_MAX, row_comm, ierr)
            IF (TRIM(op_mode) .EQ. "sigma") &
                 & op_scalar_gpu = 1.0_8 / (op_scalar_gpu * op_scalar_gpu)
         CASE ("grandsum")
            CALL gpu_ops_grandsum(ops_slot1, op_scalar_gpu)
            CALL MPI_Allreduce(MPI_IN_PLACE, op_scalar_gpu, 1, &
                 & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         END SELECT
         t1 = MPI_Wtime()
         dt_local = t1 - t0
         CALL MPI_Allreduce(dt_local, dt_op, 1, MPI_DOUBLE_PRECISION, &
              & MPI_MAX, MPI_COMM_WORLD, ierr)
         IF (timed) THEN
            op_times(NR) = dt_op
            IF (is_root) WRITE (*, '(A,1X,I0,1X,A,1X,ES12.5)') &
                 & "rep", NR, "s:", dt_op
         END IF
      END DO
      min_time = MINVAL(op_times)
      max_time = MAXVAL(op_times)
      gpu_mean = SUM(op_times) / REAL(op_reps, KIND = 8)
      gpu_med = ComputeMedian(op_times, op_reps)

      !! GPU-side result scalar. For the stateful ops (scale/increment/
      !! copy) it is the grand sum of the final state; for the reduction
      !! ops it is the last rep's global value (they are stateless).
      !! NOTE: the warmup rounds also apply the stateful ops, so their
      !! final state reflects warmup+reps applications.
      IF (TRIM(op_mode) .EQ. "scale" .OR. TRIM(op_mode) .EQ. "increment" &
           & .OR. TRIM(op_mode) .EQ. "copy") THEN
         CALL gpu_ops_grandsum(ops_slot1, op_gs1)
         gs_gpu_global = 0.0_NTREAL
         CALL MPI_Allreduce(REAL(op_gs1, KIND = NTREAL), gs_gpu_global, 1, &
              & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         op_scalar_gpu = gs_gpu_global
         op_nnz_end = INT(op_nnz_out, KIND = NTLONG)
         IF (TRIM(op_mode) .NE. "increment") op_nnz_end = op_nnz_start
      ELSE
         op_nnz_end = op_nnz_start
      END IF

      !! CPU reference: the same history (warmup + reps) through NTPoly,
      !! so timings and final states match application-for-application.
      !! NTPoly reductions return GLOBAL values on every rank (their
      !! internal allreduce) — do NOT reduce them again here.
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      t0 = MPI_Wtime()
      DO II = 1, op_warmup + op_reps
         SELECT CASE (TRIM(op_mode))
         CASE ("scale")
            CALL ScaleMatrix(matA, REAL(op_alpha, KIND = NTREAL))
         CASE ("increment")
            CALL IncrementMatrix(matB_ps, matA, alpha_in = &
                 & REAL(op_alpha, KIND = NTREAL))
         CASE ("copy")
            CALL CopyMatrix(matB_ps, matA)
         CASE ("trace")
            CALL MatrixTrace(matA, op_scalar_cpu)
         CASE ("dot")
            CALL DotMatrix(matB_ps, matA, op_scalar_cpu)
         CASE ("norm")
            op_scalar_cpu = MatrixNorm(matA)
         CASE ("sigma")
            CALL MatrixSigma(matA, op_scalar_cpu)
         CASE ("grandsum")
            CALL MatrixGrandSum(matA, op_scalar_cpu)
         END SELECT
      END DO
      t1 = MPI_Wtime()
      dt_local = t1 - t0
      CALL MPI_Allreduce(dt_local, cpu_med, 1, MPI_DOUBLE_PRECISION, &
           & MPI_MAX, MPI_COMM_WORLD, ierr)
      cpu_med = cpu_med / REAL(op_warmup + op_reps, KIND = 8)
      IF (TRIM(op_mode) == 'scale' .OR. TRIM(op_mode) == 'increment' .OR. &
           TRIM(op_mode) == 'copy') CALL MatrixGrandSum(matA, op_scalar_cpu)

      !! Parity: GPU global result vs NTPoly global result. Both are
      !! already global — no further allreduce (that would multiply by
      !! the rank count; see the MatrixGrandSum note above).
      gs_cpu_global = op_scalar_cpu
      IF (ABS(gs_cpu_global) .GT. 0.0_NTREAL) THEN
         rel_err = ABS(op_scalar_gpu - gs_cpu_global) / &
              & ABS(gs_cpu_global)
      ELSE
         rel_err = ABS(op_scalar_gpu - gs_cpu_global)
      END IF

      IF (is_root) THEN
         WRITE (*, '(A)') "-----------------"
         WRITE (*, '(A,1X,ES12.5)') "gpu min_s    :", min_time
         WRITE (*, '(A,1X,ES12.5)') "gpu median_s :", gpu_med
         WRITE (*, '(A,1X,ES12.5)') "gpu mean_s   :", gpu_mean
         WRITE (*, '(A,1X,ES12.5)') "gpu max_s    :", max_time
         WRITE (*, '(A,1X,ES12.5)') "cpu median_s :", cpu_med
         IF (cpu_med .GT. 0.0_8) WRITE (*, '(A,1X,F8.1)') &
              & "speedup      :", cpu_med / gpu_med
         WRITE (*, '(A,1X,ES24.16)') "result gpu   :", op_scalar_gpu
         WRITE (*, '(A,1X,ES24.16)') "result cpu   :", gs_cpu_global
         WRITE (*, '(A,1X,ES9.2)')  "rel_err      :", rel_err
         IF (rel_err .LT. 1.0e-12_8) THEN
            WRITE (*, '(A)') "parity       : PASS (13 sig figs)"
         ELSE
            WRITE (*, '(A)') "parity       : FAIL"
         END IF
      END IF

      !! JSON summary.
      IF (is_root) THEN
         IF (TRIM(output_file) .EQ. "") THEN
            output_file = "op_bench_results.json"
         END IF
         OPEN (NEWUNIT = unit_handle, FILE = TRIM(output_file), &
              & STATUS = 'REPLACE', ACTION = 'WRITE')
         WRITE (unit_handle, '(A)') "{"
         WRITE (unit_handle, '(A)') '  "app": "yukawa_gemm_op_bench",'
         WRITE (unit_handle, '(A,1X,I0,A)') '  "mpi_ranks":', comm_size, ","
         WRITE (unit_handle, '(A)') '  "input": "'// &
              & json_escape(TRIM(input_file))//'",'
         WRITE (unit_handle, '(A)') '  "op": "'//TRIM(op_mode)//'",'
         WRITE (unit_handle, '(A,1X,I0,A)') '  "reps":', op_reps, ","
         WRITE (unit_handle, '(A)') '  "alpha": '//real_to_json(op_alpha)//','
         WRITE (unit_handle, '(A,1X,I0,A)') '  "block_rows":', m_loc, ","
         WRITE (unit_handle, '(A,1X,I0,A)') '  "block_cols":', kk, ","
         WRITE (unit_handle, '(A,1X,I0,A)') '  "block_nnz_start":', &
              & op_nnz_start, ","
         WRITE (unit_handle, '(A,1X,I0,A)') '  "block_nnz_end":', &
              & op_nnz_end, ","
         WRITE (unit_handle, '(A)') '  "gpu_s": {'
         WRITE (unit_handle, '(A)') '    "min": '//real_to_json(min_time)//','
         WRITE (unit_handle, '(A)') '    "median": '//real_to_json(gpu_med)//','
         WRITE (unit_handle, '(A)') '    "mean": '//real_to_json(gpu_mean)//','
         WRITE (unit_handle, '(A)') '    "max": '//real_to_json(max_time)
         WRITE (unit_handle, '(A)') '  },'
         WRITE (unit_handle, '(A)') '  "cpu_median_s": '// &
              & real_to_json(cpu_med)//','
         WRITE (unit_handle, '(A)') '  "result_gpu": '// &
              & real_to_json(op_scalar_gpu)//','
         WRITE (unit_handle, '(A)') '  "result_cpu": '// &
              & real_to_json(gs_cpu_global)//','
         WRITE (unit_handle, '(A)') '  "rel_err": '//real_to_json(rel_err)
         WRITE (unit_handle, '(A)') "}"
         CLOSE (unit_handle)
         WRITE (*, '(A,1X,A)') "results     :", TRIM(output_file)
      END IF

      !! Cleanup.
      CALL DestructMatrix(matA)
      IF (TRIM(op_mode) .NE. "scale") CALL DestructMatrix(matB_ps)
      NULLIFY (mgmt_row_ptr, mgmt_col_ind, mgmt_vals)
      IF (ALLOCATED(colsum)) DEALLOCATE (colsum)
      DEALLOCATE (op_times)
      CALL DestructProcessGrid
      CALL MPI_Finalize(ierr)
      STOP
   END IF

   !! ---------------------------------------------------------------------
   !! P3 (GPU_ROADMAP.md): TRS2 density-matrix purification on resident
   !! managed slots (--solver trs2).  H = input (Hamiltonian-like) matrix;
   !! ISQ = identity so the working WH = H; X_0 = (e_max*I - WH)/
   !! (e_max - e_min) with Gershgorin bounds.  Per iteration (mirroring
   !! NTPoly's TRS2): trace(X_k) -> sigma sign; X_k2 = X_k*X_k via the 2D
   !! SUMMA cuSPARSE multiply (threshold applied); X_k = 2X_k - X_k2
   !! (sigma>0) or X_k = X_k2 (sigma<0) on managed slots; energy =
   !! <X_k, WH>; use NTPoly's convergence monitor. A reference (NTPoly's
   !! own TRS2 on identical inputs) runs first.
   IF (use_solver) THEN
      IF (solver_mode /= 'trs2' .OR. is_square .OR. use_op_bench) &
           CALL SolverAbort('use --solver trs2 --mode fixed without --op')
      IF (max_iter < 1 .OR. threshold < 0 .OR. conv_diff <= 0 .OR. chunk_work_budget < 1) &
           CALL SolverAbort('invalid iteration limit, threshold, or tolerance')
      IF (global_grid%num_process_slices /= 1) CALL SolverAbort("TRS2 requires one slice")

      CALL GET_ENVIRONMENT_VARIABLE('SLURM_JOB_ID', numbuf, STATUS=log_ios)
      IF (log_ios == 0 .AND. LEN_TRIM(numbuf) > 0) THEN
         solver_log_path = 'trs2-'//TRIM(numbuf)//'.log'
      ELSE
         IF (LEN_TRIM(output_file) == 0) output_file = 'yukawa_gemm_results.json'
         solver_log_path = TRIM(output_file)//'.trs2.log'
      END IF
      IF (is_root) THEN
         OPEN(NEWUNIT=unit_handle, FILE=TRIM(solver_log_path), STATUS='UNKNOWN', &
              POSITION='APPEND', ACTION='WRITE', IOSTAT=log_ios)
         IF (log_ios /= 0) CALL SolverAbort('cannot open solver log')
         WRITE(*,*) 'solver log: ', TRIM(solver_log_path)
      END IF
      CALL SolverStage('starting: binding runtimes')

      !! Bind the OpenACC runtime (NTPoly kernels) and the CUDA runtime
      !! (managed slots, cuSPARSE) to the SAME per-rank device BEFORE any
      !! kernel runs. Doing this late races the ACC runtime's lazy init
      !! (CUDA_ERROR_INVALID_HANDLE in MultiplyBlock).
      CALL GET_ENVIRONMENT_VARIABLE("OMPI_COMM_WORLD_LOCAL_RANK", &
           & numbuf, STATUS = ops_err)
      IF (ops_err .EQ. 0) THEN
         READ (numbuf, *) lr_local
      ELSE
         lr_local = MOD(my_rank, 4)
      END IF
      t0=MPI_Wtime()
      CALL ACC_SET_DEVICE_NUM(lr_local, ACC_DEVICE_NVIDIA)
      CALL gpu_ops_init(ops_err, ops_dev)
      IF (ops_err /= 0) CALL SolverAbort('GPU runtime/device bind failed')
      gpu_bind_seconds=MPI_Wtime()-t0
      CALL MPI_Allreduce(MPI_IN_PLACE,gpu_bind_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

      CALL SolverStage("reading and permuting input")
      !! ---- Read H (optionally load-balanced) into matWH. ----
      t0 = MPI_Wtime()
      IF (use_load_balance) THEN
         CALL ReadInputMatrix(matA_orig, TRIM(input_file))
         matrix_dim = matA_orig%actual_matrix_dimension
         CALL RANDOM_SEED(SIZE = seed_size)
         ALLOCATE (seed_array(seed_size))
         DO II = 1, seed_size
            seed_array(II) = 424242 + II
         END DO
         CALL RANDOM_SEED(PUT = seed_array)
         DEALLOCATE (seed_array)
         CALL ConstructRandomPermutation(balance_perm, matrix_dim)
         CALL PermuteMatrix(matA_orig, matWH, balance_perm)
         CALL DestructPermutation(balance_perm)
         CALL DestructMatrix(matA_orig)
      ELSE
         CALL ReadInputMatrix(matWH, TRIM(input_file))
      END IF
      tstage_read = MPI_Wtime() - t0
      matrix_dim = matWH%actual_matrix_dimension
      nnz_input = GetMatrixSize(matWH)
      m_loc = matWH%local_rows
      kk = matWH%local_columns
      n_full = matWH%logical_matrix_dimension
      IF (trace_target .LE. 0.0_NTREAL) trace_target = matrix_dim / 2.0_NTREAL
      IF (trace_target > matrix_dim) CALL SolverAbort("electron count exceeds dimension")
      IF (MeasureAsymmetry(matWH) > 1.0E-12_NTREAL) &
           CALL SolverAbort("TRS2 requires a symmetric Hamiltonian")
      my_row = matWH%process_grid%my_row
      my_col = matWH%process_grid%my_column
      row_comm = matWH%process_grid%row_comm
      col_comm = matWH%process_grid%column_comm
      pc_loc = matWH%process_grid%num_process_columns
      diag_offset = matWH%start_row - matWH%start_column
      solver_setup_t0=MPI_Wtime()
      IF (use_nccl) THEN
         t0=MPI_Wtime()
         CALL InitSolverNCCL(row_comm,col_comm)
         IF (warm_nccl) CALL WarmSolverNCCL(row_comm,col_comm)
         nccl_setup_seconds=MPI_Wtime()-t0
         CALL MPI_Allreduce(MPI_IN_PLACE,nccl_setup_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
         IF(is_root) WRITE(*,'(A,ES24.16)') 'nccl_setup_seconds: ',nccl_setup_seconds
      END IF

      !! (device already bound at the top of this branch)
      CALL GershgorinBounds(matWH, e_min, e_max)
      IF (e_max <= e_min) CALL SolverAbort("degenerate spectral bounds")

      !! X_0 = (e_max*I - WH)/(e_max - e_min) via NTPoly CPU ops (once).
      CALL ConstructEmptyMatrix(matX, matrix_dim)
      CALL CopyMatrix(matWH, matX)
      CALL ScaleMatrix(matX, -1.0_NTREAL)
      CALL ConstructEmptyMatrix(ISQ_mat, matrix_dim)
      CALL FillMatrixIdentity(ISQ_mat)
      CALL IncrementMatrix(ISQ_mat, matX, alpha_in = e_max)
      CALL ScaleMatrix(matX, 1.0_NTREAL / (e_max - e_min))

      !! NTPoly local blocks are CSC; slots always hold true CSR.
      CALL UploadSolverMatrix(matWH, solver_slot_wh)
      CALL UploadSolverMatrix(matX, solver_slot_x)
      CALL MergeMatrixLocalBlocks(matX,solver_local_x)
      IF(solver_repetitions>1) CALL gpu_ops_copy(solver_slot_initial,solver_slot_x)
      solver_setup_seconds=MPI_Wtime()-solver_setup_t0
      CALL MPI_Allreduce(MPI_IN_PLACE,solver_setup_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      IF(is_root) THEN
         WRITE(*,'(A,ES24.16)') 'gpu_bind_seconds: ',gpu_bind_seconds
         WRITE(*,'(A,ES24.16)') 'gpu_numerical_setup_seconds: ',solver_setup_seconds
      END IF

      IF (is_root) THEN
         WRITE (*, '(A)') "=== YukawaGemm (TRS2 solver) ==="
         IF(use_resident_panels) WRITE(*,*) 'panels: GPU resident CSR'
         IF(use_gpu_dense) WRITE(*,*) 'multiply: direct GPU dense accumulator'
         IF(use_device_partials) WRITE(*,*) 'partial products: GPU streaming accumulator'
         IF(use_nccl) THEN
            WRITE(*,*) "communication: NCCL device panels and scalar reductions"
         ELSE
            WRITE(*,*) "communication: MPI"
         END IF
         WRITE (*, '(A,1X,A)') "input     :", TRIM(input_file)
         WRITE (*, '(A,1X,I0)') "dimension :", matrix_dim
         WRITE (*, '(A,1X,I0)') "input nnz :", nnz_input
         WRITE (*, '(A,1X,I0)') "mpi ranks :", comm_size
         WRITE (*, '(A,1X,ES12.5)') "nelec     :", trace_target
         WRITE (*, '(A,1X,ES12.5)') "e_min     :", e_min
         WRITE (*, '(A,1X,ES12.5)') "e_max     :", e_max
         WRITE (*, '(A,1X,ES10.3)') "threshold :", threshold
         WRITE (*, '(A,1X,I0)') "max_iter  :", max_iter
         WRITE (*, '(A,1X,ES9.2)')  "conv_diff :", conv_diff
         WRITE (*, '(A)') "-----------------"
      END IF

      !! ---- CPU reference: NTPoly's own TRS2 on identical inputs. ----
      IF (is_root) THEN
         WRITE (unit_handle, '(A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES12.5)') &
              & "input", matrix_dim, "nelec", trace_target, &
              & "e_min", e_min, "e_max", e_max, "thr", threshold
         CALL FLUSH(unit_handle)
      END IF
      CALL SolverStage("starting NTPoly reference")
      CALL ConstructSolverParameters(solver_params, &
           & converge_diff_in = conv_diff, threshold_in = threshold, &
           & max_iterations_in = max_iter, be_verbose_in = .TRUE.)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      t0 = MPI_Wtime()
      IF (is_root) CALL ActivateLogger()
      CALL TRS2(matWH, ISQ_mat, trace_target, K_out, &
           & energy_value_out = solver_energy, &
           & solver_parameters_in = solver_params)
      IF (is_root) CALL DeactivateLogger()
      t1 = MPI_Wtime()
      ref_solver_seconds = t1-t0
      CALL MPI_Allreduce(MPI_IN_PLACE, ref_solver_seconds, 1, MPI_DOUBLE_PRECISION, &
           MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MatrixTrace(K_out, tr_val)
      IF (is_root) THEN
         WRITE (*, '(A)') "-----------------"
         WRITE (*, '(A,1X,ES24.16)') "ntpoly_ref_s:", ref_solver_seconds
         WRITE (*, '(A,1X,ES24.16)') "cpu energy :", solver_energy
         WRITE (*, '(A,1X,ES24.16)') "cpu trace  :", tr_val
         WRITE (*, '(A)') "-----------------"
         CALL FLUSH(6)
         WRITE (unit_handle, '(A,1X,ES24.16)') "ntpoly_ref_s:", ref_solver_seconds
         WRITE (unit_handle, '(A,1X,ES24.16)') "cpu_energy:", solver_energy
         WRITE (unit_handle, '(A,1X,ES24.16)') "cpu_trace:", tr_val
         CALL FLUSH(unit_handle)
      END IF
      reference_energy = solver_energy
      reference_trace = tr_val
      reference_norm = MatrixNorm(K_out)
      CALL DestructMatrix(ISQ_mat)

      gpu_chain_seconds=0
      DO solver_rep=1,solver_repetitions
      gpu_restore_seconds=0
      IF(solver_repetitions>1) THEN
         t0=MPI_Wtime()
         CALL gpu_ops_copy(solver_slot_x,solver_slot_initial)
         IF(.NOT. use_resident_panels) CALL RefreshSolverHost()
         gpu_restore_seconds=MPI_Wtime()-t0
         CALL MPI_Allreduce(MPI_IN_PLACE,gpu_restore_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      END IF
      solver_profile_path=solver_log_path
      IF(solver_repetitions>1) THEN
         WRITE(solver_rep_string,'(I0)') solver_rep
         solver_profile_path=TRIM(solver_log_path)//'.rep'//TRIM(solver_rep_string)
         IF(is_root) WRITE(*,'(A,I0)') 'solver_repetition: ',solver_rep
      END IF
      IF (is_root) THEN
         OPEN(NEWUNIT=profile_unit,FILE=TRIM(solver_profile_path)//'.profile.csv', &
              STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
         WRITE(profile_unit,'(A)') 'iteration,trace,row_gather,row_transpose,column_gather,'// &
              'plan_pack,chunk_prepare,spgemm_wrapper,partial_store,merge,cleanup,'// &
              'product_slot,update,host_csc,energy,wall,max_chunks,max_partial_nnz'
         FLUSH(profile_unit)
      END IF
      CALL SolverStage("starting GPU iterations")
      !! ---- GPU iteration loop. ----
      ALLOCATE (iter_times(max_iter))
      ALLOCATE(phase_history(17,max_iter))
      phase_history=0
      total_iters = 0
      energy_old = 0.0_NTREAL
      gpu_converged = .FALSE.
      CALL ConstructMonitor(gpu_monitor, automatic_in=solver_params%monitor_convergence, &
           tight_cutoff_in=conv_diff)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      CALL gpu_profile_begin()
      loop_t0 = MPI_Wtime()
      DO siter = 1, max_iter
         phase = 0
         iteration_t0 = MPI_Wtime()
         t0 = MPI_Wtime()
         !! 1. trace(X_k) -> sigma sign (resident reduction).
         CALL gpu_ops_trace(solver_slot_x, diag_offset, op_scalar_gpu)
         CALL SolverScalarSum(op_scalar_gpu)
         tr_val = REAL(op_scalar_gpu, KIND = NTREAL)
         IF (trace_target - tr_val .LT. 0.0_NTREAL) THEN
            sigma_sign = -1.0_NTREAL
         ELSE
            sigma_sign = 1.0_NTREAL
         END IF
         t_red = MPI_Wtime() - t0
         phase(1) = t_red

         !! 2. Chunked multiply X_k2 = X_k * X_k. Mid-purification fill
         !! exceeds cuSPARSE's viable range (F4), so the shared dimension
         !! is split into slabs: per-slab partial SpGEMM (cusparse_gemm_
         !! ab_partial_), partials concatenated on the host, then ONE
         !! dense-scatter merge + threshold (cusparse_ksplit_merge_) —
         !! exact, since summation completes before thresholding. The
         !! chunk count adapts to X_k's fill.
         t0 = MPI_Wtime()
         IF (use_resident_panels) THEN
            cs_ab_status=gpu_resident_square(INT(solver_slot_x,C_INT),INT(solver_slot_x2,C_INT), &
                 INT(m_loc,C_INT),INT(n_full,C_INT),INT(kk,C_INT), &
                 INT(matWH%process_grid%num_process_rows,C_INT),INT(pc_loc,C_INT), &
                 REAL(threshold,C_DOUBLE),resident_phase)
            IF(cs_ab_status/=0) CALL SolverAbort('resident panel multiply failed')
            phase(2)=resident_phase(1)
            phase(3)=resident_phase(2)
            phase(4)=resident_phase(3)
            phase(7)=resident_phase(4)
            t_mult=MPI_Wtime()-t0
            t_io=0
            nchunk=1
         ELSE
         profile_t0 = MPI_Wtime()
         row_nnz_loc = SIZE(solver_local_x%values)
         n_row_blk = m_loc
         n_col_blk = n_full
         ALLOCATE (row_blk_outer(n_col_blk + 1))
         CALL PanelFixedInt(1,solver_local_x%outer_index(2),kk,row_blk_outer(2),row_comm)
         row_blk_outer(1) = 0
         offset_acc = 0
         DO JJ = 0, pc_loc - 1
            idx_base = JJ * kk
            DO II = 1, kk
               row_blk_outer(idx_base + II + 1) = &
                    & row_blk_outer(idx_base + II + 1) + offset_acc
            END DO
            offset_acc = row_blk_outer(idx_base + kk + 1)
         END DO
         ALLOCATE (row_nnz_counts(pc_loc))
         CALL PanelFixedInt(1,[row_nnz_loc],1,row_nnz_counts,row_comm)
         row_nnz_total = SUM(row_nnz_counts)
         ALLOCATE (row_blk_inner(row_nnz_total))
         ALLOCATE (row_blk_vals(row_nnz_total))
         ALLOCATE (row_displs(pc_loc))
         row_displs(1) = 0
         DO II = 2, pc_loc
            row_displs(II) = row_displs(II-1) + row_nnz_counts(II-1)
         END DO
         CALL PanelInt(1,solver_local_x%inner_index,row_nnz_loc, &
              row_blk_inner,row_nnz_counts,row_displs,row_comm)
         CALL PanelReal(1,solver_local_x%values,row_nnz_loc, &
              row_blk_vals,row_nnz_counts,row_displs,row_comm)
         phase(2) = MPI_Wtime()-profile_t0
         profile_t0 = MPI_Wtime()
         !! Transpose the row panel CSC -> CSR (A_local).
         ALLOCATE (aloc_row_ptr(m_loc + 1))
         aloc_row_ptr = 0
         DO II = 1, row_nnz_total
            aloc_row_ptr(row_blk_inner(II) + 1) = &
                 & aloc_row_ptr(row_blk_inner(II) + 1) + 1
         END DO
         aloc_row_ptr(1) = 0
         DO II = 1, m_loc
            aloc_row_ptr(II + 1) = aloc_row_ptr(II) + aloc_row_ptr(II + 1)
         END DO
         ALLOCATE (aloc_col_ind(row_nnz_total))
         ALLOCATE (aloc_vals(row_nnz_total))
         ALLOCATE (row_count(m_loc))
         row_count = 0
         DO JJ = 1, n_full
            DO II = row_blk_outer(JJ) + 1, row_blk_outer(JJ + 1)
               r_pos = row_blk_inner(II)
               c_pos = aloc_row_ptr(r_pos) + row_count(r_pos) + 1
               aloc_col_ind(c_pos) = JJ
               aloc_vals(c_pos) = row_blk_vals(II)
               row_count(r_pos) = row_count(r_pos) + 1
            END DO
         END DO
         DEALLOCATE (row_count)
         phase(3) = MPI_Wtime()-profile_t0
         profile_t0 = MPI_Wtime()
         !! Gather the actual column panel, without a symmetry alias.
         CALL GatherSolverColumn(matX,solver_local_x, col_row_blk_outer, &
              col_row_blk_inner, col_row_blk_vals)
         b_row_ptr => col_row_blk_outer
         b_col_ind => col_row_blk_inner
         b_vals => col_row_blk_vals

         phase(4) = MPI_Wtime()-profile_t0
         profile_t0 = MPI_Wtime()
         IF (use_gpu_dense) THEN
            nchunk=1
            cs_ab_status=cusparse_gemm_ab_dense(INT(m_loc,C_INT),INT(n_full,C_INT),INT(kk,C_INT), &
                 INT(SIZE(aloc_vals),C_INT),aloc_row_ptr,aloc_col_ind,aloc_vals, &
                 INT(SIZE(b_vals),C_INT),b_row_ptr,b_col_ind,b_vals, &
                 REAL(threshold,C_DOUBLE),cs_ab_nnz_c,cs_ab_row_ptr_c,cs_ab_col_ind_c,cs_ab_val_c)
            phase(7)=MPI_Wtime()-profile_t0
         ELSE
         !! At most chunk_work_budget intermediate scalar products per slab.
         !! This bounds workspace demand; a vendor failure is still fatal.
         CALL PlanChunks(n_full, aloc_col_ind, b_row_ptr, &
              chunk_work_budget, chunk_ends, chunk_ok)
         IF (.NOT. chunk_ok) CALL SolverAbort('single-k product exceeds chunk workspace budget')
         nchunk = SIZE(chunk_ends)-1
         CALL PackChunks(m_loc, aloc_row_ptr, aloc_col_ind, aloc_vals, &
              chunk_ends, ac_rp, ac_ci, ac_val, chunk_offsets)

         phase(5) = MPI_Wtime()-profile_t0
         profile_t0 = MPI_Wtime()
         !! Per-chunk partial SpGEMM (compute on device, download the
         !! partial), concatenated on the host with offset-fixed row_ptr
         !! segments.
         IF (ALLOCATED(ck_rp)) DEALLOCATE (ck_rp)
         IF (use_device_partials) THEN
            ALLOCATE (ck_rp(1)) ! Unused dummy for the legacy merge interface.
         ELSE
            ALLOCATE (ck_rp(nchunk * (m_loc + 1)))
         END IF
         ck_rp = 0
         IF (.NOT. ALLOCATED(ck_ci)) THEN
            ALLOCATE (ck_ci(65536))
            ALLOCATE (ck_vals(65536))
         END IF
         ks_total = 0
         IF (use_device_partials) THEN
            cs_ab_status = cusparse_solver_begin(INT(m_loc,C_INT),INT(kk,C_INT))
            IF (cs_ab_status /= 0) CALL SolverAbort('device accumulator initialization failed')
         END IF
         phase(6) = MPI_Wtime()-profile_t0
         DO JJ = 1, nchunk
            profile_t0 = MPI_Wtime()
            ks_lo = chunk_ends(JJ)
            ks_hi = chunk_ends(JJ+1)
            kfull_c = ks_hi - ks_lo
            base = (JJ - 1) * (m_loc + 1)
            nnz_ac = ac_rp(base + m_loc + 1)
            nnz_bc = b_row_ptr(ks_hi + 1) - b_row_ptr(ks_lo + 1)
            IF (ALLOCATED(bc_rp_c)) DEALLOCATE (bc_rp_c)
            ALLOCATE (bc_rp_c(kfull_c + 1))
            DO II = 1, kfull_c + 1
               bc_rp_c(II) = b_row_ptr(ks_lo + II) - b_row_ptr(ks_lo + 1)
            END DO
            IF (ALLOCATED(bc_ci_c)) DEALLOCATE (bc_ci_c)
            ALLOCATE (bc_ci_c(MAX(nnz_bc, 1)))
            IF (ALLOCATED(bc_val_c)) DEALLOCATE (bc_val_c)
            ALLOCATE (bc_val_c(MAX(nnz_bc, 1)))
            IF (nnz_bc .GT. 0) THEN
               bc_ci_c(1:nnz_bc) = b_col_ind(b_row_ptr(ks_lo + 1) + 1: &
                    & b_row_ptr(ks_hi + 1))
               bc_val_c(1:nnz_bc) = b_vals(b_row_ptr(ks_lo + 1) + 1: &
                    & b_row_ptr(ks_hi + 1))
            END IF
            IF (ALLOCATED(ac_ci_c)) DEALLOCATE (ac_ci_c)
            ALLOCATE (ac_ci_c(MAX(nnz_ac, 1)))
            IF (ALLOCATED(ac_val_c)) DEALLOCATE (ac_val_c)
            ALLOCATE (ac_val_c(MAX(nnz_ac, 1)))
            IF (nnz_ac > 0) THEN
               ac_ci_c(1:nnz_ac) = ac_ci(chunk_offsets(JJ)+1:chunk_offsets(JJ+1))
               ac_val_c(1:nnz_ac) = ac_val(chunk_offsets(JJ)+1:chunk_offsets(JJ+1))
            END IF
            IF (.NOT. use_device_partials) ck_rp(base+1:base+m_loc+1) = ks_total
            phase(6) = phase(6) + MPI_Wtime()-profile_t0
            IF (nnz_ac == 0 .OR. nnz_bc == 0) CYCLE
            profile_t0 = MPI_Wtime()

            ks_nnz_p_c = 0
            cs_ab_status = cusparse_gemm_ab_partial( &
                 & INT(m_loc, KIND = C_INT), INT(kfull_c, KIND = C_INT), &
                 & INT(kk, KIND = C_INT), INT(nnz_ac, KIND = C_INT), &
                 & ac_rp(base + 1:base + m_loc + 1), &
                 & ac_ci_c(1:MAX(nnz_ac, 1)), ac_val_c(1:MAX(nnz_ac, 1)), &
                 & INT(nnz_bc, KIND = C_INT), bc_rp_c, &
                 & bc_ci_c(1:MAX(nnz_bc, 1)), bc_val_c(1:MAX(nnz_bc, 1)), &
                 & ks_nnz_p_c, ks_rp_p, ks_ci_p, ks_va_p)
            phase(7) = phase(7) + MPI_Wtime()-profile_t0
            profile_t0 = MPI_Wtime()
            IF (cs_ab_status /= 0) CALL SolverAbort('partial SpGEMM failed')
            phase(17) = phase(17) + REAL(ks_nnz_p_c,8)
            !! Device mode has already consumed the raw partial into its
            !! accumulator; it returns a count but no host output pointers.
            IF (use_device_partials) CYCLE
            IF (ks_nnz_p_c == 0) CYCLE
            IF (ks_nnz_p_c < 0 .OR. ks_nnz_p_c > HUGE(ks_total)-ks_total) &
                 CALL SolverAbort('partial product index overflow')
            CALL C_F_POINTER(ks_rp_p, ks_f_rp, [INT(m_loc + 1, KIND = C_INT)])
            CALL C_F_POINTER(ks_ci_p, ks_f_ci, [INT(ks_nnz_p_c, KIND = C_INT)])
            CALL C_F_POINTER(ks_va_p, ks_f_val, [INT(ks_nnz_p_c, KIND = C_INT)])
            CALL ReservePartials(ks_total, ks_total + ks_nnz_p_c, ck_ci, ck_vals)
            ck_rp(base + 1:base + m_loc + 1) = &
                 & ks_f_rp(1:m_loc + 1) + ks_total
            ck_ci(ks_total + 1:ks_total + ks_nnz_p_c) = &
                 & ks_f_ci(1:ks_nnz_p_c)
            ck_vals(ks_total + 1:ks_total + ks_nnz_p_c) = &
                 & ks_f_val(1:ks_nnz_p_c)
            ks_total = ks_total + ks_nnz_p_c
            phase(8) = phase(8) + MPI_Wtime()-profile_t0
         END DO
         DEALLOCATE (ac_rp, ac_ci, ac_val, chunk_ends, chunk_offsets)

         !! Merge all partials (dense scatter-add on GPU) and threshold
         !! once - exact, since summation completes before thresholding.
         profile_t0 = MPI_Wtime()
         cs_ab_status = cusparse_ksplit_merge( &
              & INT(nchunk, KIND = C_INT), INT(m_loc, KIND = C_INT), &
              & INT(kk, KIND = C_INT), INT(ks_total, KIND = C_INT), &
              & ck_rp, ck_ci(1:MAX(ks_total, 1)), &
              & ck_vals(1:MAX(ks_total, 1)), &
              & REAL(threshold, KIND = C_DOUBLE), &
              & cs_ab_nnz_c, cs_ab_row_ptr_c, cs_ab_col_ind_c, cs_ab_val_c)
         phase(9) = MPI_Wtime()-profile_t0
         END IF
         profile_t0 = MPI_Wtime()
         DEALLOCATE (row_blk_outer, row_blk_inner, row_blk_vals)
         DEALLOCATE (aloc_row_ptr, aloc_col_ind, aloc_vals)
         DEALLOCATE (row_nnz_counts, row_displs)
         DEALLOCATE (col_row_blk_outer, col_row_blk_inner, col_row_blk_vals)
         IF (cs_ab_status /= 0) CALL SolverAbort("partial-product merge failed")
         phase(10) = MPI_Wtime()-profile_t0
         t_mult = MPI_Wtime() - t0

         !! 3. Load X_k2 into its slot (host buffers -> managed).
         t0 = MPI_Wtime()
         IF (cs_ab_nnz_c .GT. 0) THEN
            CALL gpu_ops_alloc(solver_slot_x2, m_loc, cs_ab_nnz_c, &
                 & op_rp_ptr, op_ci_ptr, op_val_ptr, ops_err)
            IF (ops_err /= 0) CALL SolverAbort("product slot allocation failed")
            CALL C_F_POINTER(op_rp_ptr, mgmt_row_ptr, [m_loc + 1])
            CALL C_F_POINTER(op_ci_ptr, mgmt_col_ind, [cs_ab_nnz_c])
            CALL C_F_POINTER(op_val_ptr, mgmt_vals, [cs_ab_nnz_c])
            CALL C_F_POINTER(cs_ab_row_ptr_c, cs_ab_f_row_ptr, [m_loc + 1])
            CALL C_F_POINTER(cs_ab_col_ind_c, cs_ab_f_col_ind, [cs_ab_nnz_c])
            CALL C_F_POINTER(cs_ab_val_c, cs_ab_f_val, [cs_ab_nnz_c])
            mgmt_row_ptr = cs_ab_f_row_ptr
            mgmt_col_ind = cs_ab_f_col_ind
            mgmt_vals = cs_ab_f_val
         ELSE
            CALL gpu_ops_alloc(solver_slot_x2, m_loc, 0, op_rp_ptr, &
                 & op_ci_ptr, op_val_ptr, ops_err)
            IF (ops_err /= 0) CALL SolverAbort("empty slot allocation failed")
         END IF
         t_io = MPI_Wtime() - t0
         phase(11) = t_io

         END IF

         !! 4. Update X_k on slots: 2X - X2 (sigma>0) or X2 (sigma<0).
         t0 = MPI_Wtime()
         IF (sigma_sign .GT. 0.0_NTREAL) THEN
            CALL gpu_ops_scale(solver_slot_x, 2.0_C_DOUBLE)
            op_nnz_out = 0
            alpha_tmp = -1.0_C_DOUBLE
            CALL gpu_ops_increment_threshold(solver_slot_x, solver_slot_x2, &
                 & alpha_tmp, threshold, op_nnz_out)
            IF (op_nnz_out < 0) CALL SolverAbort("increment failed")
         ELSE
            CALL gpu_ops_copy(solver_slot_x, solver_slot_x2)
         END IF
         t_upd = MPI_Wtime() - t0
         phase(12) = t_upd

         !! 6. Energy = <X_k, WH> (resident reduction).
         t0 = MPI_Wtime()
         CALL gpu_ops_dot(solver_slot_x, solver_slot_wh, op_scalar_gpu)
         CALL SolverScalarSum(op_scalar_gpu)
         energy_new = REAL(op_scalar_gpu, KIND = NTREAL)
         IF (.NOT. IEEE_IS_FINITE(energy_new)) CALL SolverAbort("nonfinite solver energy")
         phase(14) = MPI_Wtime()-t0
         t_red = t_red + phase(14)

         !! Convert the true CSR slot back to NTPoly CSC for the next gather.
         t0 = MPI_Wtime()
         IF (use_resident_panels) THEN
            !! Read cached dimensions/nnz only; never dereference the payload.
            CALL gpu_ops_get(solver_slot_x,op_rp_ptr,op_ci_ptr,op_val_ptr,op_m_out,op_nnz_out)
         ELSE
            CALL RefreshSolverHost()
         END IF
         phase(13) = MPI_Wtime()-t0
         t_io = t_io + phase(13)

         phase(15) = MPI_Wtime()-iteration_t0
         phase(16) = REAL(nchunk,8)
         phase_history(:,siter)=phase
         CALL MPI_Reduce(phase,phase_max,SIZE(phase),MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
         IF (is_root) THEN
            WRITE(profile_unit,'(I0,17(",",ES24.16))') siter,phase_max
            FLUSH(profile_unit)
         END IF

         total_iters = siter
         iter_times(siter) = t_mult + t_red + t_upd + t_io
         IF (is_root) THEN
            WRITE (*, '(A,1X,I0,6(1X,A,1X,ES12.5))') &
                 & "it", siter, "tot:", iter_times(siter), "mult:", t_mult, &
                 & "red:", t_red, "upd:", t_upd, "io:", t_io, "E:", energy_new
            WRITE (*, '(A,1X,I0)') "xnnz:", op_nnz_out
            CALL FLUSH(6)
            WRITE (unit_handle, '(A,1X,I0,1X,ES24.16,1X,I0)') &
                 & "it", siter, energy_new, op_nnz_out
            CALL FLUSH(unit_handle)
         END IF

         !! 7. Converged?
         CALL AppendValue(gpu_monitor, energy_new-energy_old)
         gpu_converged = CheckConverged(gpu_monitor, .FALSE.)
         IF (gpu_converged) EXIT
         energy_old = energy_new
      END DO

      gpu_solver_seconds = MPI_Wtime()-loop_t0
      CALL gpu_profile_end()
      CALL MPI_Allreduce(MPI_IN_PLACE,gpu_solver_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      IF (is_root) CLOSE(profile_unit)
      !! Preserve coherent per-rank phase vectors for counterfactual timing
      !! models. Gather and file output happen after the timed/profiled loop.
      IF (is_root) THEN
         ALLOCATE(all_phase_history(17,max_iter,comm_size))
      ELSE
         ALLOCATE(all_phase_history(1,1,1))
      END IF
      CALL MPI_Gather(phase_history,17*max_iter,MPI_DOUBLE_PRECISION, &
           all_phase_history,17*max_iter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      IF (is_root) THEN
         OPEN(NEWUNIT=profile_unit,FILE=TRIM(solver_profile_path)//'.ranks.csv',STATUS='REPLACE',ACTION='WRITE')
         WRITE(profile_unit,'(A)') 'rank,iteration,trace,row_gather,row_transpose,column_gather,'// &
              'plan_pack,chunk_prepare,spgemm_wrapper,partial_store,merge,cleanup,'// &
              'product_slot,update,host_csc,energy,wall,max_chunks,max_partial_nnz'
         DO profile_rank=1,comm_size
            DO profile_iteration=1,total_iters
               WRITE(profile_unit,'(I0,",",I0,17(",",ES24.16))') profile_rank-1, &
                    profile_iteration,all_phase_history(:,profile_iteration,profile_rank)
            END DO
         END DO
         CLOSE(profile_unit)
      END IF
      DEALLOCATE(phase_history,all_phase_history)

      CALL gpu_ops_trace(solver_slot_x, diag_offset, op_scalar_gpu)
      CALL SolverScalarSum(op_scalar_gpu)
      tr_val=op_scalar_gpu
      IF(use_resident_panels) CALL RefreshSolverHost()
      parity_energy = ABS(energy_new-reference_energy)/MAX(1.0_NTREAL, ABS(reference_energy))
      parity_trace = ABS(tr_val-reference_trace)/MAX(1.0_NTREAL, ABS(reference_trace))
      IF (.NOT. IEEE_IS_FINITE(tr_val) .OR. .NOT. IEEE_IS_FINITE(reference_energy) .OR. &
           .NOT. IEEE_IS_FINITE(reference_trace) .OR. .NOT. IEEE_IS_FINITE(reference_norm)) &
           CALL SolverAbort("nonfinite parity values")
      !! Full matrix comparison catches transposition and cancellation errors
      !! that trace/energy alone can conceal. Both matrices share the same LB.
      CALL IncrementMatrix(K_out, matX, alpha_in=-1.0_NTREAL)
      parity_matrix = MatrixNorm(matX)
      parity_matrix = parity_matrix/MAX(1.0_NTREAL, reference_norm)
      IF (.NOT. IEEE_IS_FINITE(parity_matrix)) CALL SolverAbort('nonfinite matrix difference')
      IF (is_root) THEN
         WRITE (*, '(A)') "-----------------"
         WRITE (*, '(A,1X,I0)') "iterations :", total_iters
         WRITE (*, '(A,1X,ES24.16)') "gpu energy :", energy_new
         WRITE (*, '(A,1X,ES24.16)') "gpu trace  :", tr_val
         WRITE (*,'(A,1X,ES24.16)') 'gpu_solver_seconds:',gpu_solver_seconds
         IF(use_nccl) WRITE(*,'(A,ES24.16)') 'nccl_loop_plus_setup_seconds: ', &
              gpu_solver_seconds+nccl_setup_seconds
         WRITE (*,*) 'converged:', gpu_converged
         WRITE (*,*) 'relative energy/trace/matrix error:', parity_energy, parity_trace, parity_matrix
         CALL FLUSH(6)
         WRITE (unit_handle, '(A,1X,I0)') "iterations:", total_iters
         WRITE (unit_handle, '(A,1X,ES24.16)') "gpu_energy:", energy_new
         WRITE (unit_handle, '(A,1X,ES24.16)') "gpu_trace:", tr_val
         WRITE(unit_handle,*) "relative energy/trace/matrix error:", parity_energy, parity_trace, parity_matrix
         WRITE(unit_handle,*) "converged:", gpu_converged
         FLUSH(unit_handle)
      END IF

      IF (parity_energy > 1.0E-12_NTREAL .OR. parity_trace > 1.0E-12_NTREAL .OR. &
           parity_matrix > 1.0E-12_NTREAL) &
           CALL SolverAbort('CPU/GPU final matrix or scalar parity failed')
      gpu_chain_seconds=gpu_chain_seconds+gpu_restore_seconds+gpu_solver_seconds
      IF(is_root) WRITE(*,'(A,ES24.16)') 'gpu_reset_seconds: ',gpu_restore_seconds
      CALL DestructMonitor(gpu_monitor)
      DEALLOCATE(iter_times)
      END DO
      IF(is_root) THEN
         WRITE(*,'(A,ES24.16)') 'gpu_chain_seconds: ',gpu_chain_seconds
         WRITE(*,'(A,ES24.16)') 'gpu_chain_with_setup_seconds: ', &
              gpu_bind_seconds+solver_setup_seconds+gpu_chain_seconds
         CLOSE(unit_handle)
      END IF
      IF(use_resident_panels) THEN
         cs_ab_status=gpu_resident_free()
         IF(cs_ab_status/=0) CALL SolverAbort("resident cleanup failed")
      END IF
      CALL FinishSolverNCCL()
      CALL DestructMatrix(K_out)
      CALL DestructSolverParameters(solver_params)
      !! ---- Cleanup. ----
      CALL DestructMatrix(matWH)
      CALL DestructMatrix(matX)
      CALL LocalDestruct(solver_local_x)
      CALL DestructProcessGrid
      CALL MPI_Finalize(ierr)
      STOP
   END IF

   IF (use_cusparse_mpi .AND. .NOT. is_square) THEN
      !! Read the matrix using NTPoly's distributed I/O.
      t0 = MPI_Wtime()
      IF (use_load_balance) THEN
         !! Load balancing: read into a temporary, then replace the matrix
         !! with P*A*P^T where P is a fixed-seed random permutation.
         !! Rows/columns are shuffled symmetrically, so each process block
         !! gets a homogeneous mix of the global row-nnz distribution.
         !! Trace/grand-sum/nnz fingerprints are invariant (P is orthogonal).
         CALL ReadInputMatrix(matA_orig, TRIM(input_file))
         matrix_dim = matA_orig%actual_matrix_dimension
         !! Fixed seed for reproducibility (only root's shuffle matters --
         !! ConstructRandomPermutation broadcasts the result to all ranks).
         CALL RANDOM_SEED(SIZE = seed_size)
         ALLOCATE (seed_array(seed_size))
         DO II = 1, seed_size
            seed_array(II) = 424242 + II
         END DO
         CALL RANDOM_SEED(PUT = seed_array)
         DEALLOCATE (seed_array)
         CALL ConstructRandomPermutation(balance_perm, matrix_dim)
         CALL PermuteMatrix(matA_orig, matA, balance_perm)
         CALL DestructPermutation(balance_perm)
         CALL DestructMatrix(matA_orig)
         IF (is_root) WRITE (*, '(A)') "load-balance: P*A*P^T applied (random perm)"
      ELSE
         CALL ReadInputMatrix(matA, TRIM(input_file))
      END IF
      tstage_read = MPI_Wtime() - t0
      matrix_dim = matA%actual_matrix_dimension
      nnz_input = GetMatrixSize(matA)

      !! Extract local block dimensions.
      !! NTPoly stores local data in CSC format (despite Matrix_lsr name):
      !!   outer_index = col_ptr (0-based, kk+1 entries)
      !!   inner_index = row_ind (1-based LOCAL, 1..m_loc)
      !!   values      = values
      n_full = matrix_dim
      m_loc = matA%local_rows      !! N/pr (rows of local block)
      kk = matA%local_columns      !! N/pc (cols of local block)
      start_row = matA%start_row   !! first global row (0-based)
      end_row = matA%end_row       !! one past last global row

      !! Process grid info.
      pr_loc = matA%process_grid%num_process_rows
      pc_loc = matA%process_grid%num_process_columns
      ps_loc = matA%process_grid%num_process_slices
      my_row = matA%process_grid%my_row
      my_col = matA%process_grid%my_column
      row_comm = matA%process_grid%row_comm
      col_comm = matA%process_grid%column_comm
      within_comm = matA%process_grid%within_slice_comm

      !! Partner rank for symmetry exploitation (transpose position in grid).
      partner_rank = my_col * pc_loc + my_row

      !! Require square process grid (pr = pc) for the CSC<->CSR aliasing trick.
      IF (pr_loc .NE. pc_loc) THEN
         IF (is_root) WRITE (*,*) &
              & "Error: --cusparse-2d requires square process grid (pr = pc)"
         CALL DestructMatrix(matA)
         CALL DestructProcessGrid
         CALL MPI_Finalize(ierr)
         STOP
      END IF

      !! Slice requirements: plain 2D SUMMA needs ps=1; the 2.5D k-split
      !! needs ps>=2 (inner dimension split across slices) and is
      !! incompatible with --symmetric.
      IF (use_ksplit) THEN
         IF (ps_loc .LT. 2) THEN
            IF (is_root) WRITE (*,*) &
                 & "Error: --ksplit requires pslices >= 2"
            CALL DestructMatrix(matA)
            CALL DestructProcessGrid
            CALL MPI_Finalize(ierr)
            STOP
         END IF
         IF (use_symmetric) THEN
            IF (is_root) WRITE (*,*) &
                 & "Error: --ksplit cannot be combined with --symmetric"
            CALL DestructMatrix(matA)
            CALL DestructProcessGrid
            CALL MPI_Finalize(ierr)
            STOP
         END IF
      ELSE IF (ps_loc .NE. 1) THEN
         IF (is_root) WRITE (*,*) &
              & "Error: --cusparse-2d requires pslices=1 (use --ksplit)"
         CALL DestructMatrix(matA)
         CALL DestructProcessGrid
         CALL MPI_Finalize(ierr)
         STOP
      END IF

      !! Merge NTPoly's local sub-blocks into one CSC structure covering
      !! the full local block (m_loc x kk).  With pslices=1 and OMP
      !! threads=1 this is a single (1x1) block; with pslices>1 (or
      !! OMP>1) the local block is split into a (nbr x nbc) array of CSC
      !! sub-blocks, and reading only sub-block (1,1) silently drops data.
      CALL MergeLocalBlocksToCSC(matA, loc_outer, loc_inner, loc_vals, &
           & row_nnz_loc)
      t0 = MPI_Wtime()

      !! ---- Step 1: AllGather row-block along row_comm ----
      !! All ranks in the same row have the same rows but different columns.
      !! After AllGather: each rank has the full row-block (m_loc x N) in CSC.
      n_row_blk = m_loc   !! rows of row-block
      n_col_blk = n_full  !! cols of row-block

      !! AllGather outer_index (col_ptr): kk entries per rank (skip first = 0).
      ALLOCATE(row_blk_outer(n_col_blk + 1))
      CALL MPI_AllGather(loc_outer(2), kk, &
           & MPI_INTEGER, row_blk_outer(2), kk, MPI_INTEGER, &
           & row_comm, ierr)
      row_blk_outer(1) = 0
      !! Fix up offsets: each block's col_ptr is relative (starts from 0).
      !! After AllGather, add running offset to make global col_ptr.
      offset_acc = 0
      DO II = 0, pc_loc - 1
         idx_base = II * kk   !! 0-based start of block II's data
         DO JJ = 1, kk
            row_blk_outer(idx_base + JJ + 1) = &
                 & row_blk_outer(idx_base + JJ + 1) + offset_acc
         END DO
         offset_acc = row_blk_outer(idx_base + kk + 1)
      END DO

      !! AllGather nnz counts, then AllGatherv inner_index and values.
      ALLOCATE(row_nnz_counts(pc_loc))
      CALL MPI_AllGather(row_nnz_loc, 1, MPI_INTEGER, &
           & row_nnz_counts, 1, MPI_INTEGER, row_comm, ierr)
      row_nnz_total = SUM(row_nnz_counts)
      ALLOCATE(row_blk_inner(row_nnz_total))
      ALLOCATE(row_blk_vals(row_nnz_total))
      ALLOCATE(row_displs(pc_loc))
      row_displs(1) = 0
      DO II = 2, pc_loc
         row_displs(II) = row_displs(II-1) + row_nnz_counts(II-1)
      END DO
      CALL MPI_AllGatherv(loc_inner, row_nnz_loc, &
           & MPI_INTEGER, row_blk_inner, row_nnz_counts, row_displs, &
           & MPI_INTEGER, row_comm, ierr)
      CALL MPI_AllGatherv(loc_vals, row_nnz_loc, &
           & MPI_DOUBLE_PRECISION, row_blk_vals, row_nnz_counts, row_displs, &
           & MPI_DOUBLE_PRECISION, row_comm, ierr)
      DEALLOCATE (loc_outer, loc_inner, loc_vals)
      tstage_gather = MPI_Wtime() - t0

      !! row_blk is now the row-block in CSC: m_loc x N.
      !! row_blk_outer = col_ptr (0-based), row_blk_inner = row_ind (1-based
      !! LOCAL to this row, 1..m_loc), row_blk_vals = values.

      !! ---- Step 3: Transpose A_rows from CSC -> CSR -> A_local ----
      !! CSC: col_ptr (row_blk_outer, N+1, 0-based)
      !!      row_ind (row_blk_inner, row_nnz_total, 1..m_loc 1-based LOCAL)
      !!      values  (row_blk_vals)
      !! CSR: row_ptr (aloc_row_ptr, m_loc+1, 0-based)
      !!      col_ind (aloc_col_ind, row_nnz_total, 1..N 1-based GLOBAL)
      !!      values  (aloc_vals)
      t0 = MPI_Wtime()
      ALLOCATE(aloc_row_ptr(m_loc + 1))
      aloc_row_ptr = 0
      DO II = 1, row_nnz_total
         r_pos = row_blk_inner(II)   !! 1-based LOCAL row (1..m_loc)
         aloc_row_ptr(r_pos + 1) = aloc_row_ptr(r_pos + 1) + 1
      END DO
      !! Prefix sum -> 0-based row_ptr
      aloc_row_ptr(1) = 0
      DO II = 1, m_loc
         aloc_row_ptr(II + 1) = aloc_row_ptr(II) + aloc_row_ptr(II + 1)
      END DO
      !! Scatter col_ind and values
      ALLOCATE(aloc_col_ind(row_nnz_total))
      ALLOCATE(aloc_vals(row_nnz_total))
      ALLOCATE(row_count(m_loc))
      row_count = 0
      DO JJ = 1, n_full
         DO II = row_blk_outer(JJ) + 1, row_blk_outer(JJ + 1)
            r_pos = row_blk_inner(II)
            c_pos = aloc_row_ptr(r_pos) + row_count(r_pos) + 1
            aloc_col_ind(c_pos) = JJ         !! GLOBAL column (1-based)
            aloc_vals(c_pos) = row_blk_vals(II)
            row_count(r_pos) = row_count(r_pos) + 1
         END DO
      END DO
      DEALLOCATE(row_count)
      tstage_transpose = MPI_Wtime() - t0
      t0 = MPI_Wtime()

      !! ---- Step 4: Broadcast B (row-block from row my_col) along col_comm ----
      !! For symmetric A: B = A[:, my_col] in CSR (N x m_loc)
      !!               = A[my_col, :] in CSC (m_loc x N)  [transpose alias]
      !! All ranks in column my_col need the same B: the row-block from row
      !! my_col, held by rank my_col in col_comm.  Use MPI_Bcast (collective,
      !! all ranks participate -- avoids deadlock vs pairwise Sendrecv).
      !! Root (my_row == my_col): send from row_blk arrays, alias B to them.
      !! Non-root: receive into col_row_blk arrays, alias B to them.
      IF (my_row .EQ. my_col) THEN
         col_row_nnz_total = row_nnz_total
      END IF
      !! Bcast 1: nnz count
      CALL MPI_Bcast(col_row_nnz_total, 1, MPI_INTEGER, my_col, &
           & col_comm, ierr)
      !! Allocate receive buffers on non-root ranks
      IF (my_row .NE. my_col) THEN
         ALLOCATE(col_row_blk_outer(n_full + 1))
         ALLOCATE(col_row_blk_inner(col_row_nnz_total))
         ALLOCATE(col_row_blk_vals(col_row_nnz_total))
      END IF
      !! Bcast 2: col_ptr (N+1 entries, fixed size)
      IF (my_row .EQ. my_col) THEN
         CALL MPI_Bcast(row_blk_outer, n_full + 1, MPI_INTEGER, my_col, &
              & col_comm, ierr)
      ELSE
         CALL MPI_Bcast(col_row_blk_outer, n_full + 1, MPI_INTEGER, my_col, &
              & col_comm, ierr)
      END IF
      !! Bcast 3: inner_index and values (variable size)
      IF (my_row .EQ. my_col) THEN
         CALL MPI_Bcast(row_blk_inner, row_nnz_total, MPI_INTEGER, &
              & my_col, col_comm, ierr)
         CALL MPI_Bcast(row_blk_vals, row_nnz_total, MPI_DOUBLE_PRECISION, &
              & my_col, col_comm, ierr)
      ELSE
         CALL MPI_Bcast(col_row_blk_inner, col_row_nnz_total, MPI_INTEGER, &
              & my_col, col_comm, ierr)
         CALL MPI_Bcast(col_row_blk_vals, col_row_nnz_total, &
              & MPI_DOUBLE_PRECISION, my_col, col_comm, ierr)
      END IF
      !! Alias B pointers
      IF (my_row .EQ. my_col) THEN
         b_row_ptr => row_blk_outer
         b_col_ind => row_blk_inner
         b_vals    => row_blk_vals
      ELSE
         b_row_ptr => col_row_blk_outer
         b_col_ind => col_row_blk_inner
         b_vals    => col_row_blk_vals
      END IF
      tstage_bcast = MPI_Wtime() - t0
      WRITE (*, '(A,I0,A,4(A,ES12.5))') "SPROF rank ", my_rank, ":", &
           & " read_lb ", tstage_read, " gather ", tstage_gather, &
           & " transpose ", tstage_transpose, " bcast ", tstage_bcast

      !! ---- 2.5D k-split setup: slice the shared dimension across slices.
      !! Slice s works on the slab K_s = columns/rows [s*kfull+1, (s+1)*kfull]
      !! (1-based global).  A_slab = A_i·(:, K_s) as CSR (m_loc × kfull);
      !! B_slab = B_·j(K_s, :) as CSR (kfull × kk).  Partials C^(s) are
      !! summed across between_slice_comm (ranks with identical (i,j)).
      IF (use_ksplit) THEN
         my_slice = matA%process_grid%my_slice
         btw_comm = matA%process_grid%between_slice_comm
         btw_rank = matA%process_grid%between_slice_rank
         btw_size = ps_loc
         IF (MOD(n_full, ps_loc) .NE. 0) THEN
            IF (is_root) WRITE (*,*) &
                 & "Error: --ksplit requires N divisible by pslices"
            CALL DestructMatrix(matA)
            CALL DestructProcessGrid
            CALL MPI_Finalize(ierr)
            STOP
         END IF
         IF (INT(m_loc, KIND=8) * INT(kk, KIND=8) * 8_8 &
              & .GT. 64_8 * 1024_8**3) THEN
            IF (is_root) WRITE (*,*) &
                 & "Error: --ksplit dense merge buffer exceeds 64 GiB", &
                 & m_loc, "x", kk
            CALL DestructMatrix(matA)
            CALL DestructProcessGrid
            CALL MPI_Finalize(ierr)
            STOP
         END IF
         kfull = n_full / ps_loc
         ks_lo = my_slice * kfull          !! 0-based first global col of slab
         !! A_slab: keep aloc entries whose global column (1-based) lies in
         !! (ks_lo, ks_lo+kfull]; re-base columns to 1..kfull.
         nnz_a_s = 0
         DO II = 1, m_loc
            DO JJ = aloc_row_ptr(II) + 1, aloc_row_ptr(II + 1)
               IF (aloc_col_ind(JJ) .GT. ks_lo .AND. &
                    & aloc_col_ind(JJ) .LE. ks_lo + kfull) THEN
                  nnz_a_s = nnz_a_s + 1
               END IF
            END DO
         END DO
         ALLOCATE(aloc_s_row_ptr(m_loc + 1))
         ALLOCATE(aloc_s_col_ind(MAX(nnz_a_s, 1)))
         ALLOCATE(aloc_s_vals(MAX(nnz_a_s, 1)))
         aloc_s_row_ptr(1) = 0
         ks_pos = 0
         DO II = 1, m_loc
            DO JJ = aloc_row_ptr(II) + 1, aloc_row_ptr(II + 1)
               IF (aloc_col_ind(JJ) .GT. ks_lo .AND. &
                    & aloc_col_ind(JJ) .LE. ks_lo + kfull) THEN
                  ks_pos = ks_pos + 1
                  aloc_s_col_ind(ks_pos) = aloc_col_ind(JJ) - ks_lo
                  aloc_s_vals(ks_pos) = aloc_vals(JJ)
               END IF
            END DO
            aloc_s_row_ptr(II + 1) = ks_pos
         END DO
         !! B_slab: rows ks_lo+1..ks_lo+kfull of B (CSR over N rows);
         !! row_ptr re-based, col_ind stays 1-based local (1..kk).
         ks_b0 = b_row_ptr(ks_lo + 1)
         ks_b1 = b_row_ptr(ks_lo + kfull + 1)
         nnz_b_s = ks_b1 - ks_b0
         ALLOCATE(b_s_row_ptr(kfull + 1))
         ALLOCATE(b_s_col_ind(MAX(nnz_b_s, 1)))
         ALLOCATE(b_s_vals(MAX(nnz_b_s, 1)))
         DO II = 1, kfull + 1
            b_s_row_ptr(II) = b_row_ptr(ks_lo + II) - ks_b0
         END DO
         IF (nnz_b_s .GT. 0) THEN
            b_s_col_ind(1:nnz_b_s) = b_col_ind(ks_b0 + 1:ks_b1)
            b_s_vals(1:nnz_b_s) = b_vals(ks_b0 + 1:ks_b1)
         END IF
         ALLOCATE(ks_dummy_val(1))
         IF (is_root) WRITE (*, '(A,1X,I0,1X,A,1X,I0)') &
              & "ksplit    : slices", ps_loc, "k/slice", kfull
      END IF

      !! Print parameters.
      IF (is_root) THEN
         WRITE (*, '(A)') "=== YukawaGemm ==="
         WRITE (*, '(A,1X,A)') "input     :", TRIM(input_file)
         WRITE (*, '(A,1X,I0)') "dimension :", matrix_dim
         WRITE (*, '(A,1X,I0)') "input nnz :", nnz_input
         WRITE (*, '(A,1X,A)') "mode      :", TRIM(mode_string)
         WRITE (*, '(A,1X,ES10.3)') "threshold :", threshold
         WRITE (*, '(A,1X,I0)') "warmup    :", num_warmup
         WRITE (*, '(A,1X,I0)') "rounds    :", num_rounds
         WRITE (*, '(A,1X,I0)') "mpi ranks :", comm_size
         WRITE (*, '(A,1X,I0)') "grid      :", pr_loc
         WRITE (*, '(A,1X,I0)') "rows/rank :", m_loc
         WRITE (*, '(A,1X,I0)') "cols/rank :", kk
         IF (use_ksplit) THEN
            WRITE (*, '(A)') "backend   : cuSPARSE 2.5D k-split SUMMA"
         ELSE
            WRITE (*, '(A)') "backend   : cuSPARSE 2D SUMMA"
         END IF
         WRITE (*, '(A)') "-----------------"
      END IF

      ALLOCATE (gemm_times(num_rounds), rot_times(num_rounds))
      ALLOCATE (round_nnz(num_rounds))

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      !! ---- Per-rank diagnostics: A_local and B structure ----
      BLOCK
         INTEGER :: dr, max_r_nnz, min_r_nnz, max_c_nnz, min_c_nnz
         INTEGER(KIND=8) :: flop_count, max_row_flop
         INTEGER, ALLOCATABLE :: rcount(:), colcount(:), brow_nnz(:)
         INTEGER, ALLOCATABLE :: row_flop(:)
         ALLOCATE(rcount(m_loc), colcount(n_full), brow_nnz(n_full))
         !! A_local per-row nnz from row_ptr (m_loc+1, 0-based)
         DO dr = 1, m_loc
            rcount(dr) = aloc_row_ptr(dr+1) - aloc_row_ptr(dr)
         END DO
         max_r_nnz = MAXVAL(rcount)
         min_r_nnz = MINVAL(rcount)
         !! A_local per-column nnz (col_ind is 1-based GLOBAL 1..N)
         colcount = 0
         DO dr = 1, row_nnz_total
            colcount(aloc_col_ind(dr)) = colcount(aloc_col_ind(dr)) + 1
         END DO
         !! B per-row nnz from b_row_ptr (N+1, 0-based)
         DO dr = 1, n_full
            brow_nnz(dr) = b_row_ptr(dr+1) - b_row_ptr(dr)
         END DO
         max_c_nnz = MAXVAL(brow_nnz)
         min_c_nnz = HUGE(0)
         DO dr = 1, n_full
            IF (brow_nnz(dr) > 0) min_c_nnz = MIN(min_c_nnz, brow_nnz(dr))
         END DO
         !! Flop count: sum_k A_col_nnz[k] * B_row_nnz[k]
         flop_count = 0
         DO dr = 1, n_full
            flop_count = flop_count + INT(colcount(dr), KIND=8) * INT(brow_nnz(dr), KIND=8)
         END DO
         !! Max per-row flop: for each row i, sum over k in A_local[i,:]
         !!   of B_row_nnz[k]
         ALLOCATE(row_flop(m_loc))
         row_flop = 0
         DO dr = 1, m_loc
            DO JJ = aloc_row_ptr(dr)+1, aloc_row_ptr(dr+1)
               row_flop(dr) = row_flop(dr) + brow_nnz(aloc_col_ind(JJ))
            END DO
         END DO
         max_row_flop = MAXVAL(row_flop)
         WRITE (*,'(A,I0,A,I0,A,I0)') "RANK ", my_rank, &
              & " (", my_row, ",", my_col, ")"
         WRITE (*,'(A,I0,A,I0,A,I0)') "  A_local: ", m_loc, "x", n_full, &
              & " nnz=", row_nnz_total
         WRITE (*,'(A,I0,A,I0)') "  A_local row nnz: max=", max_r_nnz, &
              & " min=", min_r_nnz
         WRITE (*,'(A,I0,A,I0,A,I0)') "  B: ", n_full, "x", kk, &
              & " nnz=", col_row_nnz_total
         WRITE (*,'(A,I0,A,I0)') "  B row nnz: max=", max_c_nnz, &
              & " min=", min_c_nnz
         WRITE (*,'(A,I0,A,I15)') "  flop_count=", flop_count
         WRITE (*,'(A,I0,A,I15)') "  max_row_flop=", max_row_flop
         DEALLOCATE(rcount, colcount, brow_nnz, row_flop)
      END BLOCK

      !! Multiply loop: warmup rounds first, then timed rounds.
      !! For fixed mode (A doesn't change), gather + transpose done ONCE.
      !! Each round: just cuSPARSE compute + threshold.
      DO II = 1, num_warmup + num_rounds
         timed = II .GT. num_warmup
         NR = II - num_warmup

         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         t0 = MPI_Wtime()
         IF (use_symmetric .AND. my_row > my_col) THEN
            !! Lower triangle: receive C^T from partner (skip SpGEMM)
            tpre = MPI_Wtime()
            CALL MPI_Recv(recv_nnz, 1, MPI_INTEGER, partner_rank, 100, &
                 & MPI_COMM_WORLD, mpi_status_buf, ierr)
            cs_ab_nnz_c = INT(recv_nnz, KIND=C_INT)
            IF (.NOT. ALLOCATED(recv_row_ptr)) ALLOCATE(recv_row_ptr(m_loc + 1))
            IF (recv_nnz > 0) THEN
               IF (.NOT. ALLOCATED(recv_col_ind) .OR. &
                    & SIZE(recv_col_ind) < recv_nnz) THEN
                  IF (ALLOCATED(recv_col_ind)) DEALLOCATE(recv_col_ind)
                  ALLOCATE(recv_col_ind(recv_nnz))
               END IF
               IF (.NOT. ALLOCATED(recv_vals) .OR. &
                    & SIZE(recv_vals) < recv_nnz) THEN
                  IF (ALLOCATED(recv_vals)) DEALLOCATE(recv_vals)
                  ALLOCATE(recv_vals(recv_nnz))
               END IF
               CALL MPI_Recv(recv_row_ptr, m_loc + 1, MPI_INTEGER, &
                    & partner_rank, 101, MPI_COMM_WORLD, mpi_status_buf, ierr)
               CALL MPI_Recv(recv_col_ind, recv_nnz, MPI_INTEGER, &
                    & partner_rank, 102, MPI_COMM_WORLD, mpi_status_buf, ierr)
               CALL MPI_Recv(recv_vals, recv_nnz, MPI_DOUBLE_PRECISION, &
                    & partner_rank, 103, MPI_COMM_WORLD, mpi_status_buf, ierr)
            END IF
            cs_ab_status = 0
            trecv = MPI_Wtime() - tpre
            IF (timed .AND. NR .EQ. 1) WRITE (*, '(A,I0,A,I0,A,I0,A,ES12.5)') &
                 & "LBPROF rank ", my_rank, " (", my_row, ",", my_col, &
                 & ") lower-recv: ", trecv
         ELSE
            !! Upper triangle or diagonal: compute SpGEMM
            tpre = MPI_Wtime()
            IF (use_ksplit) THEN
               !! 2.5D k-split round: partial product over this slice's
               !! k-slab (kept in device + host persistent buffers), then
               !! exchange partials across slices and merge+threshold on GPU.
               cs_ab_nnz_a = INT(nnz_a_s, KIND=C_INT)
               cs_ab_nnz_b = INT(nnz_b_s, KIND=C_INT)
               cs_ab_status = cusparse_gemm_ab_partial( &
                    & INT(m_loc, KIND=C_INT), &
                    & INT(kfull, KIND=C_INT), &
                    & INT(kk, KIND=C_INT), &
                    & cs_ab_nnz_a, aloc_s_row_ptr, aloc_s_col_ind, aloc_s_vals, &
                    & cs_ab_nnz_b, b_s_row_ptr, b_s_col_ind, b_s_vals, &
                    & ks_nnz_p_c, ks_rp_p, ks_ci_p, ks_va_p)
               ks_nnz_p = INT(ks_nnz_p_c)
               IF (cs_ab_status .NE. 0) ks_nnz_p = 0
               tpart = MPI_Wtime() - tpre

               !! --- Exchange partials across slices (same (i,j)) ---
               tpre = MPI_Wtime()
               IF (.NOT. ALLOCATED(ks_counts)) THEN
                  ALLOCATE(ks_counts(btw_size), ks_displs(btw_size))
                  ALLOCATE(ks_rp_counts(btw_size), ks_rp_displs(btw_size))
                  ALLOCATE(ksp_row_ptr(btw_size * (m_loc + 1)))
                  ALLOCATE(ks_my_rp(m_loc + 1))
                  ks_rp_counts = m_loc + 1
                  DO JJ = 1, btw_size
                     ks_rp_displs(JJ) = (JJ - 1) * (m_loc + 1)
                  END DO
               END IF
               CALL MPI_Allgather(ks_nnz_p, 1, MPI_INTEGER, ks_counts, 1, &
                    & MPI_INTEGER, btw_comm, ierr)
               ks_total = SUM(ks_counts)
               ks_displs(1) = 0
               DO JJ = 2, btw_size
                  ks_displs(JJ) = ks_displs(JJ-1) + ks_counts(JJ-1)
               END DO
               IF (.NOT. ALLOCATED(ksp_col_ind) .OR. &
                    & SIZE(ksp_col_ind) < ks_total) THEN
                  IF (ALLOCATED(ksp_col_ind)) DEALLOCATE(ksp_col_ind)
                  IF (ALLOCATED(ksp_vals)) DEALLOCATE(ksp_vals)
                  ALLOCATE(ksp_col_ind(MAX(ks_total, 1)))
                  ALLOCATE(ksp_vals(MAX(ks_total, 1)))
               END IF
               !! Stage my row_ptr segment (zeros if the partial is empty).
               IF (ks_nnz_p .GT. 0) THEN
                  CALL C_F_POINTER(ks_rp_p, ks_f_rp, &
                       & [INT(m_loc + 1, KIND=C_INT)])
                  ks_my_rp(1:m_loc+1) = ks_f_rp(1:m_loc+1)
               ELSE
                  ks_my_rp = 0
               END IF
               CALL MPI_Allgatherv(ks_my_rp, m_loc + 1, MPI_INTEGER, &
                    & ksp_row_ptr, ks_rp_counts, ks_rp_displs, &
                    & MPI_INTEGER, btw_comm, ierr)
               IF (ks_nnz_p .GT. 0) THEN
                  CALL C_F_POINTER(ks_ci_p, ks_f_ci, &
                       & [INT(ks_nnz_p, KIND=C_INT)])
                  CALL C_F_POINTER(ks_va_p, ks_f_val, &
                       & [INT(ks_nnz_p, KIND=C_INT)])
                  CALL MPI_Allgatherv(ks_f_ci, ks_nnz_p, MPI_INTEGER, &
                       & ksp_col_ind, ks_counts, ks_displs, &
                       & MPI_INTEGER, btw_comm, ierr)
                  CALL MPI_Allgatherv(ks_f_val, ks_nnz_p, &
                       & MPI_DOUBLE_PRECISION, ksp_vals, ks_counts, &
                       & ks_displs, MPI_DOUBLE_PRECISION, btw_comm, ierr)
               ELSE
                  CALL MPI_Allgatherv(ks_my_rp, 0, MPI_INTEGER, &
                       & ksp_col_ind, ks_counts, ks_displs, &
                       & MPI_INTEGER, btw_comm, ierr)
                  CALL MPI_Allgatherv(ks_dummy_val, 0, &
                       & MPI_DOUBLE_PRECISION, ksp_vals, ks_counts, &
                       & ks_displs, MPI_DOUBLE_PRECISION, btw_comm, ierr)
               END IF
               !! Offset-fix each slice's row_ptr segment so col_ind/val
               !! entries index the concatenated arrays globally.
               DO JJ = 2, btw_size
                  ks_base = ks_displs(JJ)
                  IF (ks_base .EQ. 0) CYCLE
                  ksp_row_ptr(ks_rp_displs(JJ)+1 : ks_rp_displs(JJ)+m_loc+1) = &
                       & ksp_row_ptr(ks_rp_displs(JJ)+1 : &
                       & ks_rp_displs(JJ)+m_loc+1) + ks_base
               END DO
               texch = MPI_Wtime() - tpre

               !! --- Merge partials on GPU: scatter-add into a dense
               !! m_loc x kk buffer, threshold, compact back to CSR.
               tpre = MPI_Wtime()
               cs_ab_status = cusparse_ksplit_merge( &
                    & INT(btw_size, KIND=C_INT), &
                    & INT(m_loc, KIND=C_INT), &
                    & INT(kk, KIND=C_INT), &
                    & INT(ks_total, KIND=C_INT), &
                    & ksp_row_ptr, ksp_col_ind, ksp_vals, &
                    & REAL(threshold, KIND=C_DOUBLE), &
                    & cs_ab_nnz_c, cs_ab_row_ptr_c, cs_ab_col_ind_c, &
                    & cs_ab_val_c)
               tmerge = MPI_Wtime() - tpre
               IF (timed .AND. NR .EQ. 1) WRITE (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,4(A,ES12.5))') &
                    & "KPROF rank ", my_rank, " (", my_row, ",", my_col, &
                    & ") s", my_slice, &
                    & " pnnz ", ks_total, &
                    & " partial ", tpart, " exchange ", texch, &
                    & " merge ", tmerge, " mnnz ", REAL(cs_ab_nnz_c, 8)
            ELSE
            cs_ab_nnz_a = INT(row_nnz_total, KIND=C_INT)
            cs_ab_nnz_b = INT(col_row_nnz_total, KIND=C_INT)
            cs_ab_status = cusparse_gemm_ab_threshold( &
                 & INT(n_row_blk, KIND=C_INT), &
                 & INT(n_col_blk, KIND=C_INT), &
                 & INT(kk, KIND=C_INT), &
                 & cs_ab_nnz_a, aloc_row_ptr, aloc_col_ind, aloc_vals, &
                 & cs_ab_nnz_b, b_row_ptr, b_col_ind, b_vals, &
                 & REAL(threshold, KIND=C_DOUBLE), &
                 & cs_ab_nnz_c, cs_ab_row_ptr_c, cs_ab_col_ind_c, cs_ab_val_c)
            tspgemm = MPI_Wtime() - tpre

            IF (use_symmetric .AND. my_row < my_col) THEN
               !! Upper triangle: always send C^T to partner (even if zero)
               IF (cs_ab_status .EQ. 0 .AND. cs_ab_nnz_c > 0) THEN
                  !! Transpose C and send C^T to partner
                  CALL C_F_POINTER(cs_ab_row_ptr_c, cs_ab_f_row_ptr, [m_loc + 1])
                  CALL C_F_POINTER(cs_ab_col_ind_c, cs_ab_f_col_ind, [cs_ab_nnz_c])
                  CALL C_F_POINTER(cs_ab_val_c, cs_ab_f_val, [cs_ab_nnz_c])
                  ct_nnz = INT(cs_ab_nnz_c)
                  IF (.NOT. ALLOCATED(ct_row_ptr)) ALLOCATE(ct_row_ptr(m_loc + 1))
                  IF (.NOT. ALLOCATED(ct_col_ind) .OR. &
                       & SIZE(ct_col_ind) < ct_nnz) THEN
                     IF (ALLOCATED(ct_col_ind)) DEALLOCATE(ct_col_ind)
                     ALLOCATE(ct_col_ind(ct_nnz))
                  END IF
                  IF (.NOT. ALLOCATED(ct_vals) .OR. &
                       & SIZE(ct_vals) < ct_nnz) THEN
                     IF (ALLOCATED(ct_vals)) DEALLOCATE(ct_vals)
                     ALLOCATE(ct_vals(ct_nnz))
                  END IF
                  !! Transpose CSR C -> CSR C^T (row_ptr 0-based, col_ind 1-based)
                  ALLOCATE(ct_counts(m_loc), ct_offset(m_loc))
                  ct_counts = 0
                  DO JJ = 1, ct_nnz
                     ct_counts(INT(cs_ab_f_col_ind(JJ))) = &
                          & ct_counts(INT(cs_ab_f_col_ind(JJ))) + 1
                  END DO
                  ct_row_ptr(1) = 0_C_INT
                  DO JJ = 1, m_loc
                     ct_row_ptr(JJ + 1) = ct_row_ptr(JJ) + &
                          & INT(ct_counts(JJ), C_INT)
                  END DO
                  ct_offset = 0
                  DO lr_loc = 0, m_loc - 1
                     DO kk2 = cs_ab_f_row_ptr(lr_loc + 1) + 1, &
                             & cs_ab_f_row_ptr(lr_loc + 2)
                        r_pos = INT(cs_ab_f_col_ind(kk2))
                        c_pos = INT(ct_row_ptr(r_pos)) + ct_offset(r_pos) + 1
                        ct_col_ind(c_pos) = INT(lr_loc + 1, KIND=C_INT)
                        ct_vals(c_pos) = cs_ab_f_val(kk2)
                        ct_offset(r_pos) = ct_offset(r_pos) + 1
                     END DO
                  END DO
                  DEALLOCATE(ct_counts, ct_offset)
               ELSE
                  !! SpGEMM failed or zero product
                  ct_nnz = 0
               END IF
               ttranspose = tspgemm !! placeholder, overwritten below
               !! Always send ct_nnz (partner is waiting)
               tpre = MPI_Wtime()
               CALL MPI_Send(ct_nnz, 1, MPI_INTEGER, partner_rank, 100, &
                    & MPI_COMM_WORLD, ierr)
               IF (ct_nnz > 0) THEN
                  CALL MPI_Send(ct_row_ptr, m_loc + 1, MPI_INTEGER, &
                       & partner_rank, 101, MPI_COMM_WORLD, ierr)
                  CALL MPI_Send(ct_col_ind, ct_nnz, MPI_INTEGER, &
                       & partner_rank, 102, MPI_COMM_WORLD, ierr)
                  CALL MPI_Send(ct_vals, ct_nnz, MPI_DOUBLE_PRECISION, &
                       & partner_rank, 103, MPI_COMM_WORLD, ierr)
               END IF
               tsend = MPI_Wtime() - tpre
               IF (timed .AND. NR .EQ. 1) WRITE (*, '(A,I0,A,I0,A,I0,3(A,ES12.5))') &
                    & "LBPROF rank ", my_rank, " (", my_row, ",", my_col, &
                    & ") upper: spgemm ", tspgemm, " send ", tsend, &
                    & " (ct_nnz ", REAL(ct_nnz, 8), ")"
            END IF
            IF (my_row == my_col .AND. timed .AND. NR .EQ. 1) &
                 & WRITE (*, '(A,I0,A,I0,A,I0,A,ES12.5)') &
                 & "LBPROF rank ", my_rank, " (", my_row, ",", my_col, &
                 & ") diag: spgemm ", tspgemm
            END IF
         END IF
         t1 = MPI_Wtime()
         dt_local = t1 - t0
         CALL MPI_Allreduce(dt_local, dt_gemm, 1, MPI_DOUBLE_PRECISION, &
              & MPI_MAX, MPI_COMM_WORLD, ierr)
         dt_rot = 0.0_8
         nnz_current = cs_ab_nnz_c

         IF (timed) THEN
            gemm_times(NR) = dt_gemm
            rot_times(NR) = dt_rot
            round_nnz(NR) = nnz_current
            IF (is_root) THEN
               WRITE (*, '(A,1X,I0,1X,A,1X,ES12.5,1X,A,1X,I0)') &
                    & "round", NR, "gemm_s:", dt_gemm, "nnz:", nnz_current
            END IF
         END IF
      END DO

      !! Fingerprint via MPI_AllReduce.
      !! C is m_loc x kk in CSR. col_ind is 1-based LOCAL (1..kk).
      !! Global column (0-based) = my_col * kk + (col_ind - 1).
      !! Diagonal: global_row == global_col
      !!         -> col_ind == global_row - my_col*kk + 1
      my_col_offset = my_col * kk
      IF (use_symmetric .AND. my_row > my_col) THEN
         !! Lower triangle: use received C^T data (or zeros if empty)
         IF (cs_ab_nnz_c > 0) THEN
            cs_ab_f_row_ptr => recv_row_ptr
            cs_ab_f_col_ind => recv_col_ind(1:INT(cs_ab_nnz_c))
            cs_ab_f_val => recv_vals(1:INT(cs_ab_nnz_c))
            local_trace = 0.0_NTREAL
            local_grand_sum = 0.0_NTREAL
            DO lr_loc = 0, m_loc - 1
               global_row = start_row + lr_loc
               DO kk2 = cs_ab_f_row_ptr(lr_loc + 1) + 1, &
                       & cs_ab_f_row_ptr(lr_loc + 2)
                  local_grand_sum = local_grand_sum + cs_ab_f_val(kk2)
                  IF (cs_ab_f_col_ind(kk2) .EQ. &
                       & (global_row - my_col_offset)) THEN
                     local_trace = local_trace + cs_ab_f_val(kk2)
                  END IF
               END DO
            END DO
            local_nnz = cs_ab_nnz_c
         ELSE
            local_trace = 0.0_NTREAL
            local_grand_sum = 0.0_NTREAL
            local_nnz = 0
         END IF
      ELSE IF (cs_ab_status .EQ. 0) THEN
         CALL C_F_POINTER(cs_ab_row_ptr_c, cs_ab_f_row_ptr, [m_loc + 1])
         CALL C_F_POINTER(cs_ab_col_ind_c, cs_ab_f_col_ind, [cs_ab_nnz_c])
         CALL C_F_POINTER(cs_ab_val_c, cs_ab_f_val, [cs_ab_nnz_c])

         local_trace = 0.0_NTREAL
         local_grand_sum = 0.0_NTREAL
         DO lr_loc = 0, m_loc - 1
            global_row = start_row + lr_loc
            DO kk2 = cs_ab_f_row_ptr(lr_loc + 1) + 1, &
                    & cs_ab_f_row_ptr(lr_loc + 2)
               local_grand_sum = local_grand_sum + cs_ab_f_val(kk2)
               !! col_ind is 1-based LOCAL; diagonal when
               !! col_ind == global_row - my_col_offset
               IF (cs_ab_f_col_ind(kk2) .EQ. &
                    & (global_row - my_col_offset)) THEN
                  local_trace = local_trace + cs_ab_f_val(kk2)
               END IF
            END DO
         END DO
         local_nnz = cs_ab_nnz_c
      ELSE
         local_trace = 0.0_NTREAL
         local_grand_sum = 0.0_NTREAL
         local_nnz = 0
      END IF

      !! Fingerprint reduce over one slice only: with pslices>1 every (i,j)
      !! block exists once per slice (the k-split partials were already
      !! summed across between_slice_comm), so reducing over all ranks
      !! would count each block ps times.  within_slice_comm holds each
      !! block exactly once (and spans the whole grid when ps=1).
      CALL MPI_Allreduce(local_trace, global_trace, 1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, within_comm, ierr)
      CALL MPI_Allreduce(local_grand_sum, global_grand_sum, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, within_comm, ierr)
      CALL MPI_Allreduce(local_nnz, global_nnz, 1, MPI_INTEGER8, &
           & MPI_SUM, within_comm, ierr)

      trace_value = global_trace
      grand_sum = global_grand_sum
      asymmetry = 0.0_NTREAL
      nnz_current = global_nnz

      !! Summary statistics.
      min_time = MINVAL(gemm_times)
      max_time = MAXVAL(gemm_times)
      mean_time = SUM(gemm_times) / REAL(num_rounds, KIND = 8)
      median_time = ComputeMedian(gemm_times, num_rounds)

      IF (is_root) THEN
         WRITE (*, '(A)') "-----------------"
         WRITE (*, '(A,1X,ES12.5)') "gemm min_s   :", min_time
         WRITE (*, '(A,1X,ES12.5)') "gemm median_s:", median_time
         WRITE (*, '(A,1X,ES12.5)') "gemm mean_s  :", mean_time
         WRITE (*, '(A,1X,ES12.5)') "gemm max_s   :", max_time
         WRITE (*, '(A,1X,ES24.16)') "trace        :", trace_value
         WRITE (*, '(A,1X,ES24.16)') "grand_sum    :", grand_sum
         WRITE (*, '(A,1X,ES12.5)') "asymmetry    :", asymmetry
      END IF

      !! Write JSON summary.
      IF (is_root) THEN
         IF (TRIM(output_file) .EQ. "") THEN
            output_file = "yukawa_gemm_results.json"
         END IF
         CALL DATE_AND_TIME(date_str, time_str, zone_str)
         line = date_str(1:4)//"-"//date_str(5:6)//"-"//date_str(7:8)//"T"//&
              & time_str(1:2)//":"//time_str(3:4)//":"//time_str(5:10)//zone_str

         OPEN (NEWUNIT = unit_handle, FILE = TRIM(output_file), &
              & STATUS = 'REPLACE', ACTION = 'WRITE')
         WRITE (unit_handle, '(A)') "{"
         WRITE (unit_handle, '(A)') '  "app": "yukawa_gemm",'
         WRITE (unit_handle, '(A)') '  "timestamp": "'//TRIM(line)//'",'
         WRITE (unit_handle, '(A,1X,I0,A)') '  "mpi_ranks":', comm_size, ","
         WRITE (unit_handle, '(A)') '  "input": "'//json_escape(TRIM(input_file))//'",'
         WRITE (unit_handle, '(A,1X,I0,A)') '  "matrix_dimension":', matrix_dim, ","
         WRITE (unit_handle, '(A,1X,I0,A)') '  "input_nnz":', nnz_input, ","
         WRITE (unit_handle, '(A)') '  "mode": "'//TRIM(mode_string)//'",'
         WRITE (unit_handle, '(A)') '  "threshold": '//real_to_json(threshold)//','
         WRITE (unit_handle, '(A,1X,I0,A)') '  "warmup":', num_warmup, ","
         WRITE (unit_handle, '(A,1X,I0,A)') '  "rounds":', num_rounds, ","
      IF (use_ksplit) THEN
         backend_str = "cusparse-2.5d-ksplit"
      ELSE
         backend_str = "cusparse-2d-summa"
      END IF
      WRITE (unit_handle, '(A)') '  "backend": "'//TRIM(backend_str)//'",'
         WRITE (unit_handle, '(A)') '  "rounds_gemm_s": ['
         DO II = 1, num_rounds
            CALL FormatReal(gemm_times(II), numbuf)
            line = "    "//TRIM(numbuf)
            IF (II .LT. num_rounds) line = TRIM(line)//","
            WRITE (unit_handle, '(A)') TRIM(line)
         END DO
         WRITE (unit_handle, '(A)') '  ],'
         WRITE (unit_handle, '(A)') '  "rounds_rotate_s": ['
         DO II = 1, num_rounds
            CALL FormatReal(rot_times(II), numbuf)
            line = "    "//TRIM(numbuf)
            IF (II .LT. num_rounds) line = TRIM(line)//","
            WRITE (unit_handle, '(A)') TRIM(line)
         END DO
         WRITE (unit_handle, '(A)') '  ],'
         WRITE (unit_handle, '(A)') '  "rounds_nnz": ['
         DO II = 1, num_rounds
            WRITE (numbuf, '(I0)') round_nnz(II)
            line = "    "//TRIM(numbuf)
            IF (II .LT. num_rounds) line = TRIM(line)//","
            WRITE (unit_handle, '(A)') TRIM(line)
         END DO
         WRITE (unit_handle, '(A)') '  ],'
         WRITE (unit_handle, '(A)') '  "gemm_s": {'
         WRITE (unit_handle, '(A)') '    "min": '//real_to_json(min_time)//','
         WRITE (unit_handle, '(A)') '    "median": '//real_to_json(median_time)//','
         WRITE (unit_handle, '(A)') '    "mean": '//real_to_json(mean_time)//','
         WRITE (unit_handle, '(A)') '    "max": '//real_to_json(max_time)
         WRITE (unit_handle, '(A)') '  },'
         WRITE (unit_handle, '(A)') '  "fingerprint": {'
         WRITE (unit_handle, '(A,1X,I0,A)') '    "final_nnz":', nnz_current, ","
         WRITE (unit_handle, '(A)') '    "trace": '//real_to_json(trace_value)//','
         WRITE (unit_handle, '(A)') '    "grand_sum": '//real_to_json(grand_sum)//','
         WRITE (unit_handle, '(A)') '    "asymmetry": '//real_to_json(asymmetry)
         WRITE (unit_handle, '(A)') '  }'
         WRITE (unit_handle, '(A)') "}"
         CLOSE (unit_handle)
         WRITE (*, '(A,1X,A)') "results     :", TRIM(output_file)
      END IF

      !! Cleanup.
      CALL DestructMatrix(matA)
      DEALLOCATE (row_blk_outer, row_blk_inner, row_blk_vals)
      DEALLOCATE (aloc_row_ptr, aloc_col_ind, aloc_vals)
      IF (my_row .NE. my_col) THEN
         DEALLOCATE (col_row_blk_outer, col_row_blk_inner, col_row_blk_vals)
      END IF
      DEALLOCATE (row_nnz_counts, row_displs)
      DEALLOCATE (gemm_times, rot_times, round_nnz)
      IF (ALLOCATED(ct_row_ptr)) DEALLOCATE(ct_row_ptr)
      IF (ALLOCATED(ct_col_ind)) DEALLOCATE(ct_col_ind)
      IF (ALLOCATED(ct_vals)) DEALLOCATE(ct_vals)
      IF (ALLOCATED(recv_row_ptr)) DEALLOCATE(recv_row_ptr)
      IF (ALLOCATED(recv_col_ind)) DEALLOCATE(recv_col_ind)
      IF (ALLOCATED(recv_vals)) DEALLOCATE(recv_vals)
      IF (ALLOCATED(aloc_s_row_ptr)) DEALLOCATE(aloc_s_row_ptr)
      IF (ALLOCATED(aloc_s_col_ind)) DEALLOCATE(aloc_s_col_ind)
      IF (ALLOCATED(aloc_s_vals)) DEALLOCATE(aloc_s_vals)
      IF (ALLOCATED(b_s_row_ptr)) DEALLOCATE(b_s_row_ptr)
      IF (ALLOCATED(b_s_col_ind)) DEALLOCATE(b_s_col_ind)
      IF (ALLOCATED(b_s_vals)) DEALLOCATE(b_s_vals)
      IF (ALLOCATED(ksp_row_ptr)) DEALLOCATE(ksp_row_ptr)
      IF (ALLOCATED(ksp_col_ind)) DEALLOCATE(ksp_col_ind)
      IF (ALLOCATED(ksp_vals)) DEALLOCATE(ksp_vals)
      IF (ALLOCATED(ks_my_rp)) DEALLOCATE(ks_my_rp)
      IF (ALLOCATED(ks_dummy_val)) DEALLOCATE(ks_dummy_val)
      IF (ALLOCATED(ks_counts)) DEALLOCATE(ks_counts)
      IF (ALLOCATED(ks_displs)) DEALLOCATE(ks_displs)
      IF (ALLOCATED(ks_rp_counts)) DEALLOCATE(ks_rp_counts)
      IF (ALLOCATED(ks_rp_displs)) DEALLOCATE(ks_rp_displs)
      CALL DestructProcessGrid
      CALL MPI_Finalize(ierr)
      STOP
   END IF

  !! Read the input matrix.
  CALL ReadInputMatrix(matA, TRIM(input_file))
  matrix_dim = matA%actual_matrix_dimension
  nnz_input = GetMatrixSize(matA)

  !! Setup working matrices and the memory pool (reused across rounds, as
  !! the purification solvers do).
  CALL CopyMatrix(matA, matB)
  CALL ConstructEmptyMatrix(matC, matrix_dim)
  CALL ConstructMatrixMemoryPool(memory_pool, matB)

  is_square = TRIM(mode_string) .NE. "fixed"

  !! Direct path: bypass the distributed wrapper at np=1.
  !! The distributed MatrixMultiply_ps does gather/compose/transpose on
  !! CPU every round, thrashing managed-memory pages.  At np=1 this is pure
  !! overhead -- call the local GEMM (GemmMatrix_lsr) directly instead.
  IF (use_direct .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
     CALL LocalTranspose(matA%local_data_r(1,1), matAT_local)
     CALL LocalConstructEmpty(matC_local, matrix_dim, matrix_dim)
     sparsity_est = 4.0_NTREAL * REAL(nnz_input, KIND=NTREAL) &
          & / REAL(INT(matrix_dim, KIND=NTLONG) * matrix_dim, KIND=NTREAL)
     IF (sparsity_est .GT. 1.0_NTREAL) sparsity_est = 1.0_NTREAL
     CALL LocalConstructPool(local_pool, matrix_dim, matrix_dim, sparsity_est)
  END IF

  !! cuSPARSE path: also needs matC_local for the result, but not the
  !! transpose or memory pool (cuSPARSE manages its own workspace).
  IF (use_cusparse .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
     CALL LocalConstructEmpty(matC_local, matrix_dim, matrix_dim)
  END IF

  !! Print out the parameters.
  IF (is_root) THEN
     WRITE (*, '(A)') "=== YukawaGemm ==="
     WRITE (*, '(A,1X,A)') "input     :", TRIM(input_file)
     WRITE (*, '(A,1X,I0)') "dimension :", matrix_dim
     WRITE (*, '(A,1X,I0)') "input nnz :", nnz_input
     WRITE (*, '(A,1X,A)') "mode      :", TRIM(mode_string)
     WRITE (*, '(A,1X,ES10.3)') "threshold :", threshold
     WRITE (*, '(A,1X,I0)') "warmup    :", num_warmup
     WRITE (*, '(A,1X,I0)') "rounds    :", num_rounds
     WRITE (*, '(A,1X,I0)') "mpi ranks :", comm_size
     IF (use_cusparse .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
        WRITE (*, '(A)') "backend   : cuSPARSE (hybrid + threshold)"
     ELSE IF (use_direct .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
        WRITE (*, '(A)') "backend   : OpenACC (local GEMM, direct path)"
     ELSE
        WRITE (*, '(A)') "backend   : distributed"
     END IF
     WRITE (*, '(A)') "-----------------"
  END IF

  ALLOCATE (gemm_times(num_rounds), rot_times(num_rounds))
  ALLOCATE (round_nnz(num_rounds))

  !! Multiply loop: warmup rounds first (untimed), then timed rounds.
  DO II = 1, num_warmup + num_rounds
     timed = II .GT. num_warmup
     NR = II - num_warmup

     IF (is_square) THEN
        !! C := B*B, then rotate so B holds the newest result.
        CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
        t0 = MPI_Wtime()
        CALL MatrixMultiply(matB, matB, matC, THRESHOLD_IN = threshold, &
             & MEMORY_POOL_IN = memory_pool)
        t1 = MPI_Wtime()
        dt_local = t1 - t0
        CALL MPI_Allreduce(dt_local, dt_gemm, 1, MPI_DOUBLE_PRECISION, &
             & MPI_MAX, MPI_COMM_WORLD, ierr)

        t0 = MPI_Wtime()
        CALL CopyMatrix(matC, matB)
        t1 = MPI_Wtime()
        dt_local = t1 - t0
        CALL MPI_Allreduce(dt_local, dt_rot, 1, MPI_DOUBLE_PRECISION, &
             & MPI_MAX, MPI_COMM_WORLD, ierr)
        nnz_current = GetMatrixSize(matB)
     ELSE
        !! C := A*A, a steady state multiply.
        CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
        t0 = MPI_Wtime()
        IF (use_cusparse .AND. comm_size .EQ. 1) THEN
           !! cuSPARSE hybrid: compute full product with cuSPARSE, then
           !! threshold.  Pass matA's CSC arrays as CSR to cuSPARSE;
           !! this computes C^T in CSR = C in CSC (for symmetric C = A*A).
           cs_nnz_a = SIZE(matA%local_data_r(1,1)%values)
           cs_status = cusparse_gemm_threshold( &
                & INT(matrix_dim, KIND=C_INT), &
                & cs_nnz_a, &
                & matA%local_data_r(1,1)%outer_index, &
                & matA%local_data_r(1,1)%inner_index, &
                & matA%local_data_r(1,1)%values, &
                & REAL(threshold, KIND=C_DOUBLE), &
                & cs_nnz_c, cs_row_ptr_c, cs_col_ind_c, cs_val_c)
           IF (cs_status .EQ. 0) THEN
              !! Data is in persistent C buffers (cs_row_ptr_c, etc.).
              !! Skip the C->Fortran copy during timed rounds — just track nnz.
              !! The copy to matC_local is done once after the loop for the
              !! fingerprint.
              cs_ok = .TRUE.
           ELSE
              !! cuSPARSE failed (e.g., high fill) -- fall back
              IF (use_direct) THEN
                 CALL LocalGemm(matAT_local, matAT_local, matC_local, &
                      & IsATransposed_in = .TRUE., IsBTransposed_in = .TRUE., &
                      & threshold_in = threshold, &
                      & blocked_memory_pool_in = local_pool)
              ELSE
                 CALL MatrixMultiply(matA, matA, matC, THRESHOLD_IN = threshold, &
                      & MEMORY_POOL_IN = memory_pool)
                 !! Copy distributed result to matC_local for consistent
                 !! handling by the nnz/fingerprint code below.
                 CALL LocalDestruct(matC_local)
                 CALL LocalCopy(matC%local_data_r(1,1), matC_local)
              END IF
              cs_ok = .FALSE.
           END IF
        ELSE IF (use_direct .AND. comm_size .EQ. 1) THEN
           CALL LocalGemm(matAT_local, matAT_local, matC_local, &
                & IsATransposed_in = .TRUE., IsBTransposed_in = .TRUE., &
                & threshold_in = threshold, &
                & blocked_memory_pool_in = local_pool)
        ELSE
           CALL MatrixMultiply(matA, matA, matC, THRESHOLD_IN = threshold, &
                & MEMORY_POOL_IN = memory_pool)
        END IF
        t1 = MPI_Wtime()
        dt_local = t1 - t0
        CALL MPI_Allreduce(dt_local, dt_gemm, 1, MPI_DOUBLE_PRECISION, &
             & MPI_MAX, MPI_COMM_WORLD, ierr)
        dt_rot = 0.0_8
        IF (use_cusparse .AND. comm_size .EQ. 1 .AND. cs_ok) THEN
           nnz_current = cs_nnz_c
        ELSE IF ((use_direct .OR. use_cusparse) .AND. comm_size .EQ. 1) THEN
           nnz_current = SIZE(matC_local%values)
        ELSE
           nnz_current = GetMatrixSize(matC)
        END IF
     END IF

     IF (timed) THEN
        gemm_times(NR) = dt_gemm
        rot_times(NR) = dt_rot
        round_nnz(NR) = nnz_current
        IF (is_root) THEN
           WRITE (*, '(A,1X,I0,1X,A,1X,ES12.5,1X,A,1X,I0,1X,A,1X,F8.2)') &
                & "round", NR, "gemm_s:", dt_gemm, "nnz:", nnz_current, &
                & "nnz/row:", REAL(nnz_current, KIND = 8) / REAL(matrix_dim, KIND = 8)
        END IF
     END IF
  END DO

  !! For the direct/cusparse path, copy the local result back to the distributed
  !! matrix so the fingerprint functions (which expect Matrix_ps) work.
  !! For cuSPARSE: the data is in persistent C buffers — copy to matC_local
  !! once here, then to matC for the fingerprint.
  IF (use_cusparse .AND. comm_size .EQ. 1 .AND. .NOT. is_square .AND. cs_ok) THEN
     CALL LocalDestruct(matC_local)
     matC_local%rows = matrix_dim
     matC_local%columns = matrix_dim
     CALL C_F_POINTER(cs_row_ptr_c, cs_f_row_ptr, [matrix_dim + 1])
     ALLOCATE(matC_local%outer_index(matrix_dim + 1))
     matC_local%outer_index = cs_f_row_ptr
     CALL C_F_POINTER(cs_col_ind_c, cs_f_col_ind, [cs_nnz_c])
     ALLOCATE(matC_local%inner_index(cs_nnz_c))
     matC_local%inner_index = cs_f_col_ind
     CALL C_F_POINTER(cs_val_c, cs_f_val, [cs_nnz_c])
     ALLOCATE(matC_local%values(cs_nnz_c))
     matC_local%values = cs_f_val
     CALL LocalDestruct(matC%local_data_r(1,1))
     CALL LocalCopy(matC_local, matC%local_data_r(1,1))
  ELSE IF ((use_direct .OR. use_cusparse) .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
     CALL LocalDestruct(matC%local_data_r(1,1))
     CALL LocalCopy(matC_local, matC%local_data_r(1,1))
  END IF

  !! Compute a checksum fingerprint of the final matrix. This must match
  !! across processor counts and backends (CPU vs GPU).
  IF (is_square) THEN
     CALL MatrixTrace(matB, trace_value)
     CALL MatrixGrandSum(matB, grand_sum)
     asymmetry = MeasureAsymmetry(matB)
     nnz_current = GetMatrixSize(matB)
  ELSE
     CALL MatrixTrace(matC, trace_value)
     CALL MatrixGrandSum(matC, grand_sum)
     asymmetry = MeasureAsymmetry(matC)
     nnz_current = GetMatrixSize(matC)
  END IF

  !! Summary statistics over the timed rounds.
  min_time = MINVAL(gemm_times)
  max_time = MAXVAL(gemm_times)
  mean_time = SUM(gemm_times) / REAL(num_rounds, KIND = 8)
  median_time = ComputeMedian(gemm_times, num_rounds)

  IF (is_root) THEN
     WRITE (*, '(A)') "-----------------"
     WRITE (*, '(A,1X,ES12.5)') "gemm min_s   :", min_time
     WRITE (*, '(A,1X,ES12.5)') "gemm median_s:", median_time
     WRITE (*, '(A,1X,ES12.5)') "gemm mean_s  :", mean_time
     WRITE (*, '(A,1X,ES12.5)') "gemm max_s   :", max_time
     WRITE (*, '(A,1X,ES24.16)') "trace        :", trace_value
     WRITE (*, '(A,1X,ES24.16)') "grand_sum    :", grand_sum
     WRITE (*, '(A,1X,ES12.5)') "asymmetry    :", asymmetry
  END IF

  !! Write the JSON summary.
  IF (is_root) THEN
     IF (TRIM(output_file) .EQ. "") THEN
        output_file = "yukawa_gemm_results.json"
     END IF
     CALL DATE_AND_TIME(date_str, time_str, zone_str)
     line = date_str(1:4)//"-"//date_str(5:6)//"-"//date_str(7:8)//"T"//&
          & time_str(1:2)//":"//time_str(3:4)//":"//time_str(5:10)//zone_str

     OPEN (NEWUNIT = unit_handle, FILE = TRIM(output_file), &
          & STATUS = 'REPLACE', ACTION = 'WRITE')
     WRITE (unit_handle, '(A)') "{"
     WRITE (unit_handle, '(A)') '  "app": "yukawa_gemm",'
     WRITE (unit_handle, '(A)') '  "timestamp": "'//TRIM(line)//'",'
     WRITE (unit_handle, '(A,1X,I0,A)') '  "mpi_ranks":', comm_size, ","
     WRITE (unit_handle, '(A)') '  "input": "'//json_escape(TRIM(input_file))//'",'
     WRITE (unit_handle, '(A,1X,I0,A)') '  "matrix_dimension":', matrix_dim, ","
     WRITE (unit_handle, '(A,1X,I0,A)') '  "input_nnz":', nnz_input, ","
     WRITE (unit_handle, '(A)') '  "mode": "'//TRIM(mode_string)//'",'
     WRITE (unit_handle, '(A)') '  "threshold": '//real_to_json(threshold)//','
     WRITE (unit_handle, '(A,1X,I0,A)') '  "warmup":', num_warmup, ","
     WRITE (unit_handle, '(A,1X,I0,A)') '  "rounds":', num_rounds, ","
     WRITE (unit_handle, '(A)') '  "rounds_gemm_s": ['
     DO II = 1, num_rounds
        CALL FormatReal(gemm_times(II), numbuf)
        line = "    "//TRIM(numbuf)
        IF (II .LT. num_rounds) line = TRIM(line)//","
        WRITE (unit_handle, '(A)') TRIM(line)
     END DO
     WRITE (unit_handle, '(A)') '  ],'
     WRITE (unit_handle, '(A)') '  "rounds_rotate_s": ['
     DO II = 1, num_rounds
        CALL FormatReal(rot_times(II), numbuf)
        line = "    "//TRIM(numbuf)
        IF (II .LT. num_rounds) line = TRIM(line)//","
        WRITE (unit_handle, '(A)') TRIM(line)
     END DO
     WRITE (unit_handle, '(A)') '  ],'
     WRITE (unit_handle, '(A)') '  "rounds_nnz": ['
     DO II = 1, num_rounds
        WRITE (numbuf, '(I0)') round_nnz(II)
        line = "    "//TRIM(numbuf)
        IF (II .LT. num_rounds) line = TRIM(line)//","
        WRITE (unit_handle, '(A)') TRIM(line)
     END DO
     WRITE (unit_handle, '(A)') '  ],'
     WRITE (unit_handle, '(A)') '  "gemm_s": {'
     WRITE (unit_handle, '(A)') '    "min": '//real_to_json(min_time)//','
     WRITE (unit_handle, '(A)') '    "median": '//real_to_json(median_time)//','
     WRITE (unit_handle, '(A)') '    "mean": '//real_to_json(mean_time)//','
     WRITE (unit_handle, '(A)') '    "max": '//real_to_json(max_time)
     WRITE (unit_handle, '(A)') '  },'
     WRITE (unit_handle, '(A)') '  "fingerprint": {'
     WRITE (unit_handle, '(A,1X,I0,A)') '    "final_nnz":', nnz_current, ","
     WRITE (unit_handle, '(A)') '    "trace": '//real_to_json(trace_value)//','
     WRITE (unit_handle, '(A)') '    "grand_sum": '//real_to_json(grand_sum)//','
     WRITE (unit_handle, '(A)') '    "asymmetry": '//real_to_json(asymmetry)
     WRITE (unit_handle, '(A)') '  }'
     WRITE (unit_handle, '(A)') "}"
     CLOSE (unit_handle)
     WRITE (*, '(A,1X,A)') "results     :", TRIM(output_file)
  END IF

  !! Cleanup
  CALL DestructMatrix(matA)
  CALL DestructMatrix(matB)
  CALL DestructMatrix(matC)
  CALL DestructMatrixMemoryPool(memory_pool)
  IF ((use_direct .OR. use_cusparse) .AND. comm_size .EQ. 1 .AND. .NOT. is_square) THEN
     IF (use_direct) THEN
        CALL LocalDestruct(matAT_local)
        CALL LocalDestructPool(local_pool)
     END IF
     CALL LocalDestruct(matC_local)
  END IF
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)

CONTAINS
  SUBROUTINE RefreshSolverHost()
    CALL gpu_ops_get(solver_slot_x,op_rp_ptr,op_ci_ptr,op_val_ptr,op_m_out,op_nnz_out)
    CALL C_F_POINTER(op_rp_ptr, mgmt_row_ptr, [m_loc + 1])
    CALL C_F_POINTER(op_ci_ptr, mgmt_col_ind, [MAX(op_nnz_out, 1)])
    CALL C_F_POINTER(op_val_ptr, mgmt_vals, [MAX(op_nnz_out, 1)])
    CALL TransposeCSR(m_loc, kk, mgmt_row_ptr, mgmt_col_ind, mgmt_vals, &
    loc_outer, loc_inner, loc_vals)
    CALL LocalDestruct(solver_local_x)
    solver_local_x%rows = m_loc
    solver_local_x%columns = kk
    CALL MOVE_ALLOC(loc_outer, solver_local_x%outer_index)
    CALL MOVE_ALLOC(loc_inner, solver_local_x%inner_index)
    CALL MOVE_ALLOC(loc_vals, solver_local_x%values)
    CALL SplitMatrixToLocalBlocks(matX,solver_local_x)
  END SUBROUTINE RefreshSolverHost

  SUBROUTINE SolverStage(message)
    CHARACTER(*), INTENT(IN) :: message
    IF (is_root) THEN
      WRITE(unit_handle,*) 'stage: ', message
      CALL FLUSH(unit_handle)
    END IF
  END SUBROUTINE

  SUBROUTINE SolverAbort(message)
    CHARACTER(*), INTENT(IN) :: message
    INTEGER :: abort_ierr
    WRITE (*,*) 'TRS2 error on rank', my_rank, ': ', message
    CALL MPI_Abort(MPI_COMM_WORLD, 1, abort_ierr)
    ERROR STOP
  END SUBROUTINE

  SUBROUTINE UploadSolverMatrix(mat, slot)
    TYPE(Matrix_ps), INTENT(IN) :: mat
    INTEGER(C_INT), INTENT(IN) :: slot
    TYPE(Matrix_lsr) :: local
    INTEGER, ALLOCATABLE :: rp(:), ci(:)
    REAL(C_DOUBLE), ALLOCATABLE :: va(:)
    INTEGER(C_INT) :: status, nr, nz
    TYPE(C_PTR) :: pr, pi, pv
    INTEGER(C_INT), POINTER :: fr(:), fi(:)
    REAL(C_DOUBLE), POINTER :: fv(:)
    CALL MergeMatrixLocalBlocks(mat,local)
    CALL TransposeCSR(mat%local_columns, mat%local_rows, &
         local%outer_index, local%inner_index, &
         local%values, rp, ci, va)
    nr = mat%local_rows
    nz = SIZE(ci)
    CALL gpu_ops_alloc(slot, nr, nz, pr, pi, pv, status)
    IF (status /= 0) CALL SolverAbort('initial slot allocation failed')
    CALL C_F_POINTER(pr, fr, [nr+1])
    CALL C_F_POINTER(pi, fi, [MAX(nz,1)])
    CALL C_F_POINTER(pv, fv, [MAX(nz,1)])
    fr = rp
    fi(1:nz) = ci
    fv(1:nz) = va
    CALL LocalDestruct(local)
  END SUBROUTINE

  SUBROUTINE GatherSolverColumn(mat, full_local, rp, ci, va)
    TYPE(Matrix_ps), INTENT(IN) :: mat
    TYPE(Matrix_lsr), INTENT(IN) :: full_local
    INTEGER, ALLOCATABLE, INTENT(OUT) :: rp(:), ci(:)
    REAL(NTREAL), ALLOCATABLE, INTENT(OUT) :: va(:)
    INTEGER :: nr, nc, np, nz, total, rankrow, c, j, g, p, err
    INTEGER, ALLOCATABLE :: counts(:), displs(:), ptrs(:,:), inds(:), cursor(:)
    REAL(NTREAL), ALLOCATABLE :: vals(:)
    nr = mat%local_rows
    nc = mat%local_columns
    np = mat%process_grid%num_process_rows
    nz = full_local%outer_index(nc+1)
    ALLOCATE(counts(np), displs(np), ptrs(nc+1,np))
    CALL PanelFixedInt(2,[nz],1,counts,col_comm)
    displs(1) = 0
    DO j = 2, np
      displs(j) = displs(j-1)+counts(j-1)
    END DO
    total = SUM(counts)
    ALLOCATE(inds(total), vals(total), rp(n_full+1), ci(total), va(total), cursor(n_full))
    CALL PanelFixedInt(2,full_local%outer_index,nc+1,ptrs,col_comm)
    CALL PanelInt(2,full_local%inner_index,nz,inds,counts,displs,col_comm)
    CALL PanelReal(2,full_local%values,nz,vals,counts,displs,col_comm)
    rp = 0
    DO rankrow = 1, np
      DO j = displs(rankrow)+1, displs(rankrow)+counts(rankrow)
        g = (rankrow-1)*nr+inds(j)
        rp(g+1) = rp(g+1)+1
      END DO
    END DO
    DO g = 1, n_full
      rp(g+1) = rp(g+1)+rp(g)
    END DO
    cursor = rp(1:n_full)
    DO rankrow = 1, np
      DO c = 1, nc
        DO j = ptrs(c,rankrow)+1, ptrs(c+1,rankrow)
          p = displs(rankrow)+j
          g = (rankrow-1)*nr+inds(p)
          cursor(g) = cursor(g)+1
          ci(cursor(g)) = c
          va(cursor(g)) = vals(p)
        END DO
      END DO
    END DO
  END SUBROUTINE

  !> Merge NTPoly's local sub-blocks into one CSC structure covering the
  !> whole local block (local_rows x local_columns).  With pslices=1 and
  !> OMP threads=1 there is a single (1x1) block; with pslices>1 (or
  !> OMP>1) the local block is split into a (nbr x nbc) array of CSC
  !> sub-blocks (uniform tiles, remainder in the last), each with
  !> sub-local 1-based row indices.  Gathering a single sub-block would
  !> silently drop data.
  SUBROUTINE MergeLocalBlocksToCSC(mat, full_outer, full_inner, full_vals, &
       & full_nnz)
    TYPE(Matrix_ps), INTENT(IN) :: mat
    INTEGER, ALLOCATABLE, INTENT(OUT) :: full_outer(:) !! local_columns+1, 0-based col_ptr
    INTEGER, ALLOCATABLE, INTENT(OUT) :: full_inner(:) !! 1-based local rows (1..local_rows)
    REAL(NTREAL), ALLOCATABLE, INTENT(OUT) :: full_vals(:)
    INTEGER, INTENT(OUT) :: full_nnz
    !! Local
    INTEGER :: nbr, nbc, srows, scols, bs_r, bs_c, r0, c0
    INTEGER :: I, J, JJ, II2, pos, total_cols
    INTEGER, ALLOCATABLE :: counts(:), cursor(:)

    nbr = mat%process_grid%number_of_blocks_rows
    nbc = mat%process_grid%number_of_blocks_columns
    srows = mat%local_rows / nbr
    scols = mat%local_columns / nbc
    total_cols = mat%local_columns

    full_nnz = 0
    DO J = 1, nbc
       DO I = 1, nbr
          full_nnz = full_nnz + SIZE(mat%local_data_r(I,J)%values)
       END DO
    END DO

    !! Column counts across all row sub-blocks, then 0-based col_ptr.
    ALLOCATE(counts(total_cols + 1))
    counts = 0
    DO J = 1, nbc
       bs_c = scols
       IF (J .EQ. nbc) bs_c = total_cols - scols*(nbc-1)
       c0 = (J-1)*scols
       DO JJ = 1, bs_c
          DO I = 1, nbr
             counts(c0 + JJ + 1) = counts(c0 + JJ + 1) + &
                  & mat%local_data_r(I,J)%outer_index(JJ+1) - &
                  & mat%local_data_r(I,J)%outer_index(JJ)
          END DO
       END DO
    END DO
    ALLOCATE(full_outer(total_cols + 1))
    full_outer(1) = 0
    DO JJ = 1, total_cols
       full_outer(JJ+1) = full_outer(JJ) + counts(JJ+1)
    END DO
    DEALLOCATE(counts)

    ALLOCATE(full_inner(MAX(full_nnz, 1)))
    ALLOCATE(full_vals(MAX(full_nnz, 1)))
    ALLOCATE(cursor(total_cols))
    cursor = 0
    DO J = 1, nbc
       bs_c = scols
       IF (J .EQ. nbc) bs_c = total_cols - scols*(nbc-1)
       c0 = (J-1)*scols
       DO I = 1, nbr
          bs_r = srows
          IF (I .EQ. nbr) bs_r = mat%local_rows - srows*(nbr-1)
          r0 = (I-1)*srows
          DO JJ = 1, bs_c
             DO II2 = mat%local_data_r(I,J)%outer_index(JJ)+1, &
                  & mat%local_data_r(I,J)%outer_index(JJ+1)
                pos = full_outer(c0+JJ) + cursor(c0+JJ) + 1
                full_inner(pos) = mat%local_data_r(I,J)%inner_index(II2) + r0
                full_vals(pos) = mat%local_data_r(I,J)%values(II2)
                cursor(c0+JJ) = cursor(c0+JJ) + 1
             END DO
          END DO
       END DO
    END DO
    DEALLOCATE(cursor)
  END SUBROUTINE MergeLocalBlocksToCSC

  !> Read a matrix from either a matrix market or a NTPoly binary file,
  !> dispatching on the file extension.
  SUBROUTINE ReadInputMatrix(this, file_name)
    USE PSMatrixModule, ONLY : Matrix_ps, ConstructMatrixFromMatrixMarket, &
         & ConstructMatrixFromBinary
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    CHARACTER(LEN = *), INTENT(IN) :: file_name
    IF (LEN_TRIM(file_name) .GT. 4) THEN
       IF (file_name(LEN_TRIM(file_name) - 3:LEN_TRIM(file_name)) .EQ. ".bin") THEN
          CALL ConstructMatrixFromBinary(this, file_name)
          RETURN
       END IF
    END IF
    CALL ConstructMatrixFromMatrixMarket(this, file_name)
  END SUBROUTINE ReadInputMatrix

  !> Read a NTPoly binary file directly and build CSR arrays.
  !> Each rank reads the full matrix (replicated).  Since the Yukawa matrix
  !> is symmetric, CSR = CSC and we can pass these arrays directly to
  !> cuSPARSE as both A_local and B.
  SUBROUTINE ReadBinaryCSR(filename, n, nnz, row_ptr, col_ind, vals)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: n
    INTEGER(NTLONG), INTENT(OUT) :: nnz
    INTEGER, ALLOCATABLE, INTENT(OUT) :: row_ptr(:)  ! 0-based, size n+1
    INTEGER, ALLOCATABLE, INTENT(OUT) :: col_ind(:)  ! 1-based, size nnz
    REAL(NTREAL), ALLOCATABLE, INTENT(OUT) :: vals(:)  ! size nnz
    !! Local
    INTEGER :: unit, ios, i
    INTEGER :: matrix_info(3)
    INTEGER(NTLONG) :: total_values
    INTEGER, ALLOCATABLE :: coo_col(:), coo_row(:)
    REAL(NTREAL), ALLOCATABLE :: coo_val(:)
    INTEGER, ALLOCATABLE :: row_count(:)
    INTEGER :: r, pos

    OPEN(NEWUNIT=unit, FILE=TRIM(filename), ACCESS='STREAM', &
         & FORM='UNFORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
       WRITE(*,*) "Error: cannot open ", TRIM(filename)
       CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    END IF

    !! Header: 3 x int32 (rows, cols, complex_flag) + 1 x int64 (nnz)
    READ(unit) matrix_info(1), matrix_info(2), matrix_info(3)
    READ(unit) total_values
    n = matrix_info(1)
    nnz = total_values

    !! Read COO triplets (col, row, val), 1-based indices
    ALLOCATE(coo_col(INT(nnz)), coo_row(INT(nnz)), coo_val(INT(nnz)))
    READ(unit) (coo_col(i), coo_row(i), coo_val(i), i=1,INT(nnz))
    CLOSE(unit)

    !! Build CSR: row_ptr (0-based), col_ind (1-based), vals
    ALLOCATE(row_ptr(n+1), col_ind(INT(nnz)), vals(INT(nnz)), row_count(n))
    row_count = 0
    DO i = 1, INT(nnz)
       row_count(coo_row(i)) = row_count(coo_row(i)) + 1
    END DO
    row_ptr(1) = 0
    DO i = 1, n
       row_ptr(i+1) = row_ptr(i) + row_count(i)
    END DO
    !! Fill col_ind and vals (preserving row order from COO)
    row_count = 0
    DO i = 1, INT(nnz)
       r = coo_row(i)
       pos = row_ptr(r) + row_count(r) + 1
       col_ind(pos) = coo_col(i)
       vals(pos) = coo_val(i)
       row_count(r) = row_count(r) + 1
    END DO

    DEALLOCATE(coo_col, coo_row, coo_val, row_count)
  END SUBROUTINE ReadBinaryCSR

  !> Median of an array of times (sorts a copy, small arrays).
  FUNCTION ComputeMedian(values, n) RESULT(med)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: values(n)
    REAL(8) :: med
    REAL(8) :: work(n)
    INTEGER :: i, j
    REAL(8) :: v
    work = values
    DO i = 2, n
       v = work(i)
       j = i - 1
       DO WHILE (j .GE. 1)
          IF (work(j) <= v) EXIT
          work(j + 1) = work(j)
          j = j - 1
       END DO
       work(j + 1) = v
    END DO
    IF (MOD(n, 2) .EQ. 1) THEN
       med = work((n + 1) / 2)
    ELSE
       med = 0.5_8 * (work(n / 2) + work(n / 2 + 1))
    END IF
  END FUNCTION ComputeMedian

  !> Render a double as a JSON number.
  FUNCTION real_to_json(value) RESULT(out)
    REAL(8), INTENT(IN) :: value
    CHARACTER(LEN = 40) :: out
    WRITE (out, '(ES24.16E3)') value
    out = ADJUSTL(out)
  END FUNCTION real_to_json

  !> Render a double into a string for array elements.
  SUBROUTINE FormatReal(value, out)
    REAL(8), INTENT(IN) :: value
    CHARACTER(LEN = *), INTENT(OUT) :: out
    WRITE (out, '(ES24.16E3)') value
    out = ADJUSTL(out)
  END SUBROUTINE FormatReal

  !> Escape the characters that are not legal inside a JSON string.
  FUNCTION json_escape(str_in) RESULT(str_out)
    CHARACTER(LEN = *), INTENT(IN) :: str_in
    CHARACTER(LEN = LEN(str_in)) :: str_out
    INTEGER :: i, k
    str_out = ""
    k = 0
    DO i = 1, LEN_TRIM(str_in)
       k = k + 1
       IF (str_in(i:i) .EQ. '"') THEN
          str_out(k:k + 1) = '\"'
          k = k + 1
       ELSE
          str_out(k:k) = str_in(i:i)
       END IF
    END DO
  END FUNCTION json_escape

  !> Print the command line usage.
  SUBROUTINE PrintUsage()
    WRITE (*, '(A)') "Usage: mpirun -np P yukawa_gemm --input FILE [options]"
    WRITE (*, '(A)') ""
    WRITE (*, '(A)') "Required:"
    WRITE (*, '(A)') "  --input FILE      Matrix market (.mtx) or NTPoly binary (.bin) matrix."
    WRITE (*, '(A)') ""
    WRITE (*, '(A)') "Options:"
    WRITE (*, '(A)') "  --mode MODE       square (default) or fixed."
    WRITE (*, '(A)') "  --rounds N        Timed multiply rounds (default 10)."
    WRITE (*, '(A)') "  --warmup N        Untimed warmup rounds (default 3)."
    WRITE (*, '(A)') "  --threshold T     Sparsity threshold (default 1e-4)."
    WRITE (*, '(A)') "  --out FILE        JSON summary file."
    WRITE (*, '(A)') "  --prows R         Explicit process grid rows."
    WRITE (*, '(A)') "  --pcols C         Explicit process grid columns."
    WRITE (*, '(A)') "  --pslices S       Explicit process grid slices."
    WRITE (*, '(A)') "  --symmetric       Exploit A symmetry: only compute upper"
    WRITE (*, '(A)') "                    triangle, transpose for lower (2D SUMMA)."
    WRITE (*, '(A)') "  --load-balance    Shuffle rows/cols P*A*P^T (random perm, fixed"
    WRITE (*, '(A)') "                    seed) to homogenize per-block sparsity."
    WRITE (*, '(A)') "  --convert [FILE]  Convert --input to NTPoly binary format and exit."
    WRITE (*, '(A)') "  --direct          Bypass distributed wrapper (np=1, fixed mode only)."
    WRITE (*, '(A)') "  --cusparse        Use cuSPARSE SpGEMM + threshold (np=1, fixed mode only)."
    WRITE (*, '(A)') "  --cusparse-mpi    Multi-GPU cuSPARSE with 1D row decomposition (fixed mode only)."
    WRITE (*, '(A)') "  --ksplit          2.5D k-split: split the shared dimension across"
    WRITE (*, '(A)') "                    pslices, exchange partials, merge on GPU"
    WRITE (*, '(A)') "                    (fixed mode, with --cusparse-2d)."
    WRITE (*, '(A)') "  --op MODE         GPU-resident op benchmark:"
    WRITE (*, '(A)') "                    scale | increment | copy | trace |"
    WRITE (*, '(A)') "                    dot | norm | sigma | grandsum."
    WRITE (*, '(A)') "  --nccl           NCCL solver communication, single node, up to 4 GPUs."
    WRITE (*, '(A)') "  --solver-reps N  Repeat solves with GPU X0 and communicators reused."
    WRITE (*, '(A)') "  --resident-panels Keep CSR payloads, panel assembly and products on GPU."
    WRITE (*, '(A)') "  --gpu-dense       Direct GPU sparse multiply into a dense accumulator."
    WRITE (*, '(A)') "  --device-partials Keep solver partial products in the GPU accumulator."
    WRITE (*, '(A)') "  --nccl-warmup     Warm collective paths before solving; report setup separately."
    WRITE (*, '(A)') "  --solver MODE     Solver benchmark: trs2 (fixed mode)."
    WRITE (*, '(A)') "  --nelec N         Electron count for --solver (default dim/2)."
    WRITE (*, '(A)') "  --chunk-work N   Max intermediate products per solver slab (default 4000000)."
    WRITE (*, '(A)') "  --max-iter N      Max solver iterations (default 60)."
    WRITE (*, '(A)') "  --conv-diff R     Energy convergence threshold (default 1e-6)."
    WRITE (*, '(A)') "  --reps N          Repetitions for --op (default 20)."
    WRITE (*, '(A)') "  --alpha R         Scale/increment factor (default 1.01)."
  END SUBROUTINE PrintUsage
END PROGRAM YukawaGemm
