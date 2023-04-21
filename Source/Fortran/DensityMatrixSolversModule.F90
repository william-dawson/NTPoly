!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Quantum Chemistry Systems using Purification.
MODULE DensityMatrixSolversModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE FermiOperatorModule, ONLY : ComputeDenseFOE
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteListElement, WriteHeader, &
       & EnterSubLog, ExitSubLog
  USE NTMPIModule
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixMultiply, &
       & DotMatrix, MatrixTrace, ScaleMatrix, SimilarityTransform
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix, PrintMatrixInformation, FillMatrixIdentity, &
       & TransposeMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PM
  PUBLIC :: TRS2
  PUBLIC :: TRS4
  PUBLIC :: HPCP
  PUBLIC :: DenseDensity
  PUBLIC :: ScaleAndFold
  PUBLIC :: EnergyDensityMatrix
  PUBLIC :: McWeenyStep
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  !> Based on the PM algorithm presented in \cite palser1998canonical
  SUBROUTINE PM(H, ISQ, trace, K, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WH
    TYPE(Matrix_ps) :: IMat
    TYPE(Matrix_ps) :: ISQT
    TYPE(Matrix_ps) :: X_k, X_k2, X_k3, Temp
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: factor
    REAL(NTREAL) :: lambda,alpha,alpha1,alpha2
    REAL(NTREAL) :: a1,a2,a3
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: trace_value2
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II, JJ
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="PM")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("palser1998canonical")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    ALLOCATE(sigma_array(params%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(K, H)
    CALL ConstructEmptyMatrix(WH, H)
    CALL ConstructEmptyMatrix(X_k, H)
    CALL ConstructEmptyMatrix(X_k2, H)
    CALL ConstructEmptyMatrix(X_k3, H)
    CALL ConstructEmptyMatrix(Temp, H)
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in=params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WH, e_min, e_max)

    !! Initialize
    CALL CopyMatrix(WH, X_k)

    !! Compute lambda
    CALL MatrixTrace(X_k, trace_value)
    lambda = trace_value/X_k%actual_matrix_dimension

    !! Compute alpha
    alpha1 = trace/(e_max-lambda)
    alpha2 = (X_k%actual_matrix_dimension-trace)/(lambda-e_min)
    alpha = MIN(alpha1,alpha2)

    factor = -alpha/X_k%actual_matrix_dimension

    CALL ScaleMatrix(X_k, factor)
    factor = (alpha*lambda+trace)/X_k%actual_matrix_dimension
    CALL IncrementMatrix(IMat, X_k, alpha_in=factor)

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute X_k2
       CALL MatrixMultiply(X_k, X_k, X_k2, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       !! Compute X_k3
       CALL MatrixMultiply(X_k, X_k2, X_k3, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       !! Compute X_k - X_k2
       CALL CopyMatrix(X_k, Temp)
       CALL IncrementMatrix(X_k2, Temp, &
            & alpha_in=-1.0_NTREAL, threshold_in=params%threshold)

       !! Compute Sigma
       CALL MatrixTrace(Temp, trace_value)
       CALL DotMatrix(Temp, X_k, trace_value2)
       !! If we hit 0 exact convergence, avoid a division by zero.
       IF (trace_value .LE. TINY(trace_value)) THEN
          sigma_array(II) = 1.0_NTREAL
       ELSE
          sigma_array(II) = trace_value2/trace_value
       END IF

       IF (sigma_array(II) .GT. 0.5_NTREAL) THEN
          a1 = 0.0_NTREAL
          a2 = 1.0_NTREAL + 1.0_NTREAL/sigma_array(II)
          a3 = -1.0_NTREAL/sigma_array(II)
       ELSE
          a1 = (1.0_NTREAL - 2.0_NTREAL*sigma_array(II)) &
               & / (1.0_NTREAL - sigma_array(II))
          a2 = (1.0_NTREAL + sigma_array(II)) &
               & / (1.0_NTREAL - sigma_array(II))
          a3 = -1.0_NTREAL/(1.0_NTREAL - sigma_array(II))
       END IF

       !! Update X_k
       CALL ScaleMatrix(X_k, a1)
       CALL IncrementMatrix(X_k2, X_k, &
            & alpha_in=a2, threshold_in=params%threshold)
       CALL IncrementMatrix(X_k3, X_k, &
            & alpha_in=a3, threshold_in=params%threshold)

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WH, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement("Energy Value", VALUE=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = II-1
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(X_k, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(X_k3)
    CALL DestructMatrix(Temp)
    CALL DestructMatrix(IMat)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO II = 1, params%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO JJ = 1, total_iterations
             IF (sigma_array(JJ) .GT. 0.5_NTREAL) THEN
                zero_value = ((1.0_NTREAL + sigma_array(JJ)) &
                     & *zero_value**2) - (zero_value**3)
                zero_value = zero_value/sigma_array(JJ)
             ELSE
                zero_value = ((1.0_NTREAL - 2.0_NTREAL* &
                     & sigma_array(JJ))*zero_value) &
                     & + ((1.0_NTREAL + sigma_array(JJ))* &
                     & zero_value**2) - (zero_value**3)
                zero_value = zero_value/(1.0_NTREAL - sigma_array(JJ))
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value - 0.5_NTREAL) .LT. params%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = lambda - &
            & (H%actual_matrix_dimension*midpoint - trace) &
            & /alpha
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(params)
  END SUBROUTINE PM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  !> Based on the TRS2 algorithm presented in \cite niklasson2002.
  SUBROUTINE TRS2(H, ISQ, trace, K, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WH
    TYPE(Matrix_ps) :: IMat
    TYPE(Matrix_ps) :: ISQT
    TYPE(Matrix_ps) :: X_k, X_k2, Temp
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II, JJ
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="TRS2")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("niklasson2002expansion")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    ALLOCATE(sigma_array(params%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(K, H)
    CALL ConstructEmptyMatrix(WH, H)
    CALL ConstructEmptyMatrix(X_k, H)
    CALL ConstructEmptyMatrix(X_k2, H)
    CALL ConstructEmptyMatrix(Temp, H)
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in=params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WH, e_min, e_max)

    !! Initialize
    CALL CopyMatrix(WH, X_k)
    CALL ScaleMatrix(X_k, -1.0_NTREAL)
    CALL IncrementMatrix(IMat, X_k, alpha_in=e_max)
    CALL ScaleMatrix(X_k, 1.0_NTREAL/(e_max - e_min))

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute Sigma
       CALL MatrixTrace(X_k, trace_value)
       IF (trace - trace_value .LT. 0.0_NTREAL) THEN
          sigma_array(II) = -1.0_NTREAL
       ELSE
          sigma_array(II) = 1.0_NTREAL
       END IF

       !! Compute X_k2
       CALL MatrixMultiply(X_k, X_k, X_k2, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       !! Update X_k
       IF (sigma_array(II) .GT. 0.0_NTREAL) THEN
          CALL ScaleMatrix(X_k, 2.0_NTREAL)
          CALL IncrementMatrix(X_k2, X_k, &
               & alpha_in=-1.0_NTREAL, threshold_in=params%threshold)
       ELSE
          CALL CopyMatrix(X_k2,X_k)
       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WH, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement("Energy Value", VALUE=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = II - 1
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(X_k, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(Temp)
    CALL DestructMatrix(IMat)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO II = 1, params%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO JJ = 1, total_iterations
             IF (sigma_array(JJ) .LT. 0.0_NTREAL) THEN
                zero_value = zero_value*zero_value
             ELSE
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. params%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(params)
  END SUBROUTINE TRS2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  !> Based on the TRS4 algorithm presented in \cite niklasson2002
  SUBROUTINE TRS4(H, ISQ, trace, K, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN)  :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: sigma_min = 0.0_NTREAL
    REAL(NTREAL), PARAMETER :: sigma_max = 6.0_NTREAL
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WH
    TYPE(Matrix_ps) :: IMat
    TYPE(Matrix_ps) :: ISQT
    TYPE(Matrix_ps) :: X_k, X_k2, Fx_right, GX_right, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    REAL(NTREAL) :: tempfx,tempgx
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II, JJ
    INTEGER :: total_iterations
    REAL(NTREAL) :: trace_fx, trace_gx

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="TRS4")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("niklasson2002expansion")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    ALLOCATE(sigma_array(params%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(K, H)
    CALL ConstructEmptyMatrix(WH, H)
    CALL ConstructEmptyMatrix(X_k, H)
    CALL ConstructEmptyMatrix(X_k2, H)
    CALL ConstructEmptyMatrix(TempMat, H)
    CALL ConstructEmptyMatrix(Fx_right, H)
    CALL ConstructEmptyMatrix(Gx_right, H)
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in=params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WH,e_min,e_max)

    !! Initialize
    CALL CopyMatrix(WH,X_k)
    CALL ScaleMatrix(X_k, -1.0_NTREAL)
    CALL IncrementMatrix(IMat, X_k, alpha_in=e_max)
    CALL ScaleMatrix(X_k, 1.0_NTREAL/(e_max - e_min))

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute X_k2
       CALL MatrixMultiply(X_k, X_k, X_k2, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       !! Compute Fx_right
       CALL CopyMatrix(X_k2, Fx_right)
       CALL ScaleMatrix(Fx_right, -3.0_NTREAL)
       CALL IncrementMatrix(X_k, Fx_right, alpha_in=4.0_NTREAL)
       !! Compute Gx_right
       CALL CopyMatrix(IMat, Gx_right)
       CALL IncrementMatrix(X_k, Gx_right, alpha_in=-2.0_NTREAL)
       CALL IncrementMatrix(X_k2, Gx_right)

       !! Compute Traces
       CALL DotMatrix(X_k2, Fx_right, trace_fx)
       CALL DotMatrix(X_k2, Gx_right, trace_gx)

       !! Avoid Overflow
       IF (ABS(trace_gx) .LT. 1.0e-14_NTREAL) THEN
          EXIT
       END IF

       !! Compute Sigma
       sigma_array(II) = (trace - trace_fx)/trace_gx

       !! Update The Matrix
       IF (sigma_array(II) .GT. sigma_max) THEN
          CALL CopyMatrix(X_k, TempMat)
          CALL ScaleMatrix(TempMat, 2.0_NTREAL)
          CALL IncrementMatrix(X_k2, TempMat, alpha_in=-1.0_NTREAL)
       ELSE IF (sigma_array(II) .LT. sigma_min) THEN
          CALL CopyMatrix(X_k2, TempMat)
       ELSE
          CALL ScaleMatrix(Gx_right, sigma_array(II))
          CALL IncrementMatrix(Fx_right, Gx_right)
          CALL MatrixMultiply(X_k2, Gx_right, TempMat, &
               & threshold_in=params%threshold, memory_pool_in=pool)
       END IF

       CALL IncrementMatrix(TempMat, X_k, alpha_in=-1.0_NTREAL)
       CALL CopyMatrix(TempMat, X_k)

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WH, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement("Energy Value", VALUE=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = II - 1
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(X_k, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(Fx_right)
    CALL DestructMatrix(Gx_right)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(IMat)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO II = 1, params%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO JJ = 1, total_iterations
             IF (sigma_array(JJ) .GT. sigma_max) THEN
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             ELSE IF (sigma_array(JJ) .LT. sigma_min) THEN
                zero_value = zero_value*zero_value
             ELSE
                tempfx = (zero_value*zero_value) * &
                     & (4.0_NTREAL*zero_value - &
                     &  3.0_NTREAL*zero_value*zero_value)
                tempgx = (zero_value*zero_value) * (1.0_NTREAL - zero_value) &
                     & * (1.0_NTREAL - zero_value)
                zero_value = tempfx + sigma_array(JJ)*tempgx
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. params%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(params)
  END SUBROUTINE TRS4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  !> Based on the algorithm presented in \cite truflandier2016communication.
  SUBROUTINE HPCP(H, ISQ, trace, K, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WH
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: IMat
    TYPE(Matrix_ps) :: ISQT
    TYPE(Matrix_ps) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value, norm_value2
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II, JJ
    INTEGER :: total_iterations
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="HPCP")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("truflandier2016communication")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    ALLOCATE(sigma_array(params%max_iterations))

    matrix_dimension = H%actual_matrix_dimension

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(K, H)
    CALL ConstructEmptyMatrix(WH, H)
    CALL ConstructEmptyMatrix(TempMat, H)
    CALL ConstructEmptyMatrix(D1, H)
    CALL ConstructEmptyMatrix(DH, H)
    CALL ConstructEmptyMatrix(DDH, H)
    CALL ConstructEmptyMatrix(D2DH, H)
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in=params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the initial matrix.
    CALL GershgorinBounds(WH, e_min, e_max)
    CALL MatrixTrace(WH, mu)
    mu = mu/matrix_dimension
    sigma_bar = (matrix_dimension - trace)/matrix_dimension
    sigma = 1.0_NTREAL - sigma_bar
    beta = sigma/(e_max - mu)
    beta_bar = sigma_bar/(mu - e_min)
    beta_1 = sigma
    beta_2 = MIN(beta,beta_bar)

    !! Initialize
    CALL CopyMatrix(IMat, D1)
    CALL ScaleMatrix(D1, beta_1)
    CALL CopyMatrix(IMat, TempMat)
    CALL ScaleMatrix(TempMat, mu)
    CALL IncrementMatrix(WH, TempMat, -1.0_NTREAL)
    CALL ScaleMatrix(TempMat, beta_2)
    CALL IncrementMatrix(TempMat, D1)
    trace_value = 0.0_NTREAL

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF

    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    norm_value2 = norm_value
    energy_value = 0.0_NTREAL
    DO II = 1, params%max_iterations
       !! Compute the hole matrix DH
       CALL CopyMatrix(D1, DH)
       CALL IncrementMatrix(IMat, DH, alpha_in=-1.0_NTREAL)
       CALL ScaleMatrix(DH, -1.0_NTREAL)

       !! Compute DDH, as well as convergence check
       CALL MatrixMultiply(D1, DH, DDH, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL MatrixTrace(DDH, trace_value)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL MatrixMultiply(D1, DDH, D2DH, &
            & threshold_in=params%threshold, memory_pool_in=pool)

       !! Compute Sigma
       CALL MatrixTrace(D2DH, sigma_array(II))
       sigma_array(II) = sigma_array(II)/trace_value

       CALL CopyMatrix(D1, TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementMatrix(D2DH, D1, alpha_in=2.0_NTREAL)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementMatrix(DDH, D1, &
            & alpha_in=-1.0_NTREAL*2.0_NTREAL*sigma_array(II))

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(D1, WH, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement("Energy Value", VALUE=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = II - 1
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II)
       CALL PrintMatrixInformation(D1)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(D1, D1, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(D1, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(D1)
    CALL DestructMatrix(DH)
    CALL DestructMatrix(DDH)
    CALL DestructMatrix(D2DH)
    CALL DestructMatrix(IMat)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO II = 1, params%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO JJ = 1, total_iterations
             zero_value = zero_value + &
                  & 2.0_NTREAL*((zero_value**2)*(1.0_NTREAL-zero_value) &
                  & - sigma_array(JJ)* &
                  & zero_value*(1.0_NTREAL-zero_value))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. params%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = mu + (beta_1 - midpoint)/beta_2
    END IF
    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(params)
  END SUBROUTINE HPCP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the Scale and Fold
  !> method. Based on the method of \cite rubensson2011nonmonotonic .
  !> Note that for this method, you must provide the value of the homo and
  !> lumo gap. It is not necessary for these to be accurate, but give a
  !> conservative value.
  SUBROUTINE ScaleAndFold(H, ISQ, trace, K, &
       & homo, lumo, energy_value_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> A conservative estimate of the highest occupied eigenvalue.
    REAL(NTREAL), INTENT(IN) :: homo
    !> A conservative estimate of the lowest unoccupied eigenvalue.
    REAL(NTREAL), INTENT(IN) :: lumo
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Matrices
    TYPE(Matrix_ps) :: WH
    TYPE(Matrix_ps) :: IMat
    TYPE(Matrix_ps) :: ISQT
    TYPE(Matrix_ps) :: X_k, X_k2, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: Beta, BetaBar, alpha
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", VALUE="Scale and Fold")
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("rubensson2011nonmonotonic")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(K, H)
    CALL ConstructEmptyMatrix(WH, H)
    CALL ConstructEmptyMatrix(X_k, H)
    CALL ConstructEmptyMatrix(X_k2, H)
    CALL ConstructEmptyMatrix(TempMat, H)
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in=params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WH, e_min, e_max)

    !! Initialize
    CALL CopyMatrix(WH, X_k)
    CALL ScaleMatrix(X_k, -1.0_NTREAL)
    CALL IncrementMatrix(IMat, X_k, alpha_in=e_max)
    CALL ScaleMatrix(X_k, 1.0_NTREAL/(e_max - e_min))
    Beta = (e_max - lumo) / (e_max - e_min)
    BetaBar = (e_max - homo) / (e_max - e_min)

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO II = 1, params%max_iterations
       !! Determine the path
       CALL MatrixTrace(X_k, trace_value)
       IF (trace_value .GT. trace) THEN
          alpha = 2.0/(2.0 - Beta)
          CALL ScaleMatrix(X_k, alpha)
          CALL IncrementMatrix(IMat, X_k, alpha_in=(1.0_NTREAL-alpha))
          CALL MatrixMultiply(X_k, X_k, X_k2, &
               & threshold_in=params%threshold, memory_pool_in=pool)
          CALL CopyMatrix(X_k2, X_k)
          Beta = (alpha * Beta + 1 - alpha)**2
          BetaBar = (alpha * BetaBar + 1 - alpha)**2
       ELSE
          alpha = 2.0/(1.0 + BetaBar)
          CALL MatrixMultiply(X_k, X_k, X_k2, &
               & threshold_in=params%threshold, memory_pool_in=pool)
          CALL ScaleMatrix(X_k, 2*alpha)
          CALL IncrementMatrix(X_k2, X_k, alpha_in=-1.0_NTREAL*alpha**2)
          Beta = 2.0 * alpha * Beta - alpha**2 * Beta**2
          BetaBar = 2.0 * alpha * BetaBar - alpha**2 * BetaBar ** 2
       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WH, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
          CALL EnterSubLog
          CALL WriteElement("Energy Value", VALUE=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total Iterations", VALUE=II)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(X_k, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(IMat)
    CALL DestructMatrixMemoryPool(pool)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE ScaleAndFold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix using a dense routine.
  SUBROUTINE DenseDensity(H, ISQ, trace, K, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    REAL(NTREAL) :: chemical_potential, energy_value

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Call the unified routine.
    CALL ComputeDenseFOE(H, ISQ, trace, K, energy_value_out=energy_value, &
         & chemical_potential_out=chemical_potential, &
         & solver_parameters_in=params)

    !! Optional out variables.
    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF
    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = chemical_potential
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)
  END SUBROUTINE DenseDensity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the energy-weighted density matrix.
  SUBROUTINE EnergyDensityMatrix(H, D, ED, threshold_in)
    !> The matrix to compute from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The density matrix.
    TYPE(Matrix_ps), INTENT(IN) :: D
    !> The energy-weighted density matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: ED
    !> Threshold for flushing small values (default = 0).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Handling Optional Parameters
    REAL(NTREAL) :: threshold

    !! Optional Parameters
    IF (PRESENT(threshold_in)) THEN
       threshold = threshold_in
    ELSE
       threshold = 0.0_NTREAL
    END IF

    !! EDM = DM * H * DM
    CALL SimilarityTransform(H, D, D, ED, threshold_in=threshold)

  END SUBROUTINE EnergyDensityMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one McWeeny Step DOut = 3*DSD - 2*DSDSD
  SUBROUTINE McWeenyStep(D, DOut, S_in, threshold_in)
    !> The density matrix.
    TYPE(Matrix_ps), INTENT(IN) :: D
    !> The resulting purified matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: DOut
    !> The overlap matrix (optional)
    TYPE(Matrix_ps), INTENT(IN), OPTIONAL :: S_in
    !> Threshold for flushing small values (default = 0).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Handling Optional Parameters
    REAL(NTREAL) :: threshold
    !! Local Variables
    TYPE(Matrix_ps) :: DS, DSD
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(threshold_in)) THEN
       threshold = threshold_in
    ELSE
       threshold = 0.0_NTREAL
    END IF

    !! Form the matrix DS
    IF (PRESENT(S_in)) THEN
       CALL MatrixMultiply(D, S_in, DS, &
            & threshold_in=threshold, memory_pool_in=pool)
    ELSE
       CALL CopyMatrix(D, DS)
    END IF

    !! Compute
    CALL MatrixMultiply(DS, D, DSD, & 
         & threshold_in=threshold, memory_pool_in=pool)
    CALL MatrixMultiply(DS, DSD, DOut, alpha_in=-2.0_NTREAL, &
         & threshold_in=threshold, memory_pool_in=pool)
    CALL IncrementMatrix(DSD, DOut, alpha_in=3.0_NTREAL)

    !! Cleanup
    CALL DestructMatrix(DS)
    CALL DestructMatrix(DSD)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE McWeenyStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule
