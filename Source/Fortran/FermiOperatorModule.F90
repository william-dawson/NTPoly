!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Density Matrix Using the Fermi Operator Expansion
MODULE FermiOperatorModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, &
       & EnterSubLog, ExitSubLog, WriteListElement, WriteComment
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, SimilarityTransform, &
       & IncrementMatrix, MatrixNorm, ScaleMatrix, DotMatrix, MatrixTrace, &
       & MatrixDiagonalScale
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & FillMatrixFromTripletList, GetMatrixTripletList, &
       & TransposeMatrix, ConjugateMatrix, DestructMatrix, FilterMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, CopyMatrix, &
       & GatherMatrixTripletList, GetMatrixSize
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & PrintParameters, DestructSolverParameters, CopySolverParameters, &
       & ConstructSolverParameters
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & DestructTripletList, CopyTripletList
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeDenseFOE
  PUBLIC :: WOM_GC
  PUBLIC :: WOM_C
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix using a dense routine.
  SUBROUTINE ComputeDenseFOE(H, ISQ, trace, K, inv_temp_in, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The trace of the density matrix (usually the number of electrons)
    REAL(NTREAL), INTENT(IN) :: trace
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The inverse temperature for smearing in a.u. (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: inv_temp_in
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    REAL(NTREAL) :: inv_temp
    LOGICAL :: do_smearing
    !! Local Variables
    TYPE(Matrix_ps) :: ISQT, WH
    TYPE(Matrix_ps) :: WD
    TYPE(Matrix_ps) :: vecs, vecsT, vals, Temp
    TYPE(MatrixMemoryPool_p) :: pool
    TYPE(TripletList_r) :: tlist
    TYPE(TripletList_c) :: tlist_c
    REAL(NTREAL) :: chemical_potential, energy_value
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: eigs, occ
    REAL(NTREAL) :: sval, sv, occ_temp
    REAL(NTREAL) :: left, right, homo, lumo
    INTEGER :: num_eigs
    INTEGER :: II, JJ

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF
    IF (PRESENT(inv_temp_in)) THEN
       inv_temp = inv_temp_in
       do_smearing = .TRUE.
    ELSE
       do_smearing = .FALSE.
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       IF (do_smearing) THEN
          CALL WriteElement(key = "Method", VALUE = "Dense FOE")
          CALL WriteElement(key = "Inverse Temperature", VALUE = inv_temp)
       ELSE
          CALL WriteElement(key = "Method", VALUE = "Dense Step Function")
       END IF
       CALL PrintParameters(params)
    END IF

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL MatrixMultiply(ISQ, H, Temp, &
         & threshold_in = params%threshold, memory_pool_in = pool)
    CALL MatrixMultiply(Temp, ISQT, WH, &
         & threshold_in = params%threshold, memory_pool_in = pool)

    !! Perform the eigendecomposition
    CALL EigenDecomposition(WH, vals, &
         & eigenvectors_in = vecs, solver_parameters_in = params)

    !! Gather the eigenvalues on to every process
    CALL GatherMatrixTripletList(vals, tlist)

    !! Put them in an array for simplicity
    num_eigs = H%actual_matrix_dimension
    ALLOCATE(eigs(num_eigs))
    eigs = 0
    DO II = 1, tlist%CurrentSize
       eigs(II) = tlist%DATA(II)%point_value
    END DO

    !! Compute MU By Bisection
    IF (do_smearing) THEN
       ALLOCATE(occ(num_eigs))
       left = MINVAL(eigs)
       right = MAXVAL(eigs)
       DO JJ = 1, 10*params%max_iterations 
          chemical_potential = left + (right - left) / 2 
          DO II = 1, num_eigs
             sval = eigs(II) - chemical_potential
             ! occ(II) = 0.5_NTREAL * (1.0_NTREAL - ERF(inv_temp * sval))
             occ(II) = 1.0_NTREAL / (1.0_NTREAL + EXP(inv_temp * sval))
          END DO
          sv = SUM(occ)
          IF (ABS(trace - sv) .LT. 1E-8_NTREAL) THEN
             EXIT
          ELSE IF (SV > trace) THEN
             right = chemical_potential
          ELSE
             left = chemical_potential
          END IF
       END DO
    ELSE
       JJ = 1
       homo = eigs(FLOOR(trace))
       lumo = eigs(FLOOR(trace) + 1)
       occ_temp = FLOOR(TRACE) + 1 - trace
       chemical_potential = homo + occ_temp * 0.5_NTREAL * (lumo - homo)
    END IF

    !! Write out result of chemical potential search
    IF (params%be_verbose) THEN
       CALL WriteHeader("Chemical Potential Search")
       CALL EnterSubLog
       CALL WriteElement(key = "Potential", VALUE = chemical_potential)
       CALL WriteElement(key = "Iterations", VALUE = JJ)
       CALL ExitSubLog
    END IF

    !! Map - note that we store the square root of the occupation numbers
    energy_value = 0.0_NTREAL
    DO II = 1, tlist%CurrentSize
       IF (.NOT. do_smearing) THEN
          IF (tlist%DATA(II)%index_column .LE. FLOOR(trace)) THEN
             energy_value = energy_value + tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = 1.0_NTREAL
          ELSE IF (tlist%DATA(II)%index_column .EQ. CEILING(trace)) THEN
             occ_temp = trace - FLOOR(trace)
             energy_value = energy_value + &
                  & occ_temp * tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = SQRT(occ_temp)
          ELSE
             tlist%DATA(II)%point_value = 0.0_NTREAL
          ENDIF
       ELSE
          sval = tlist%DATA(II)%point_value - chemical_potential
          ! occ_temp = 0.5_NTREAL * (1.0_NTREAL - ERF(inv_temp * sval))
          occ_temp = 1.0_NTREAL / (1.0_NTREAL + EXP(inv_temp * sval))
          energy_value = energy_value + &
             & occ_temp * tlist%DATA(II)%point_value
          IF (occ_temp .LT. 0) THEN  ! for safety
             tlist%DATA(II)%point_value = 0
          ELSE
             tlist%DATA(II)%point_value = SQRT(occ_temp)
          END IF
       END IF
    END DO

    !! Scale the eigenvectors
    IF (vecs%is_complex) THEN
       CALL CopyTripletList(tlist, tlist_c)
       CALL MatrixDiagonalScale(vecs, tlist_c)
    ELSE
       CALL MatrixDiagonalScale(vecs, tlist)
    END IF
    CALL FilterMatrix(vecs, params%threshold)

    !! Multiply Back Together
    CALL TransposeMatrix(vecs, vecsT)
    CALL ConjugateMatrix(vecsT)
    CALL MatrixMultiply(vecs, vecsT, WD, &
         & threshold_in = params%threshold)

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(ISQT, WD, Temp, &
         & threshold_in = params%threshold, memory_pool_in = pool)
    CALL MatrixMultiply(Temp, ISQ, K, &
         & threshold_in = params%threshold, memory_pool_in = pool)

    !! Optional out variables.
    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF
    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = chemical_potential
    END IF

    !! Cleanup
    CALL DestructMatrix(WH)
    CALL DestructMatrix(WD)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(vecs)
    CALL DestructMatrix(vecst)
    CALL DestructMatrix(vals)
    CALL DestructMatrix(temp)
    CALL DestructTripletList(tlist)
    CALL DestructTripletList(tlist_c)
    CALL DestructMatrixMemoryPool(pool)
    IF (ALLOCATED(occ)) THEN
       DEALLOCATE(occ)
    END IF
    IF (ALLOCATED(eigs)) THEN
       DEALLOCATE(eigs)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE ComputeDenseFOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix according to the wave operator minization
  !! method at fixed chemical potential.
  SUBROUTINE WOM_GC(H, ISQ, K, chemical_potential, inv_temp, &
       & energy_value_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The inverse temperature for smearing in a.u.
    REAL(NTREAL), INTENT(IN) :: inv_temp
    !> The chemical potential.
    REAL(NTREAL), INTENT(IN) :: chemical_potential
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Print out details
    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key = "Method", VALUE = "WOM_GC")
       CALL WriteElement(key = "Inverse Temperature", VALUE = inv_temp)
       CALL WriteElement(key = "Chemical Potential", &
            & VALUE = chemical_potential)
       CALL PrintParameters(params)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       CALL WOM_Implementation(H, ISQ, K, inv_temp, params, &
            & mu_in = chemical_potential, energy_value_out = energy_value_out)
    ELSE
       CALL WOM_Implementation(H, ISQ, K, inv_temp, params, &
            & mu_in = chemical_potential)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructSolverParameters(params)
  END SUBROUTINE WOM_GC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix according to the wave operator minization
  !! method at fixed number of electrons.
  SUBROUTINE WOM_C(H, ISQ, K, trace, inv_temp, &
       & energy_value_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The inverse temperature for smearing in a.u.
    REAL(NTREAL), INTENT(IN) :: inv_temp
    !> The target trace.
    REAL(NTREAL), INTENT(IN) :: trace
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    !! Print out details
    IF (params%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key = "Method", VALUE = "WOM_C")
       CALL WriteElement(key = "Inverse Temperature", VALUE = inv_temp)
       CALL WriteElement(key = "Target Trace", VALUE = trace)
       CALL PrintParameters(params)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       CALL WOM_Implementation(H, ISQ, K, inv_temp, params, &
            & trace_in = trace, energy_value_out = energy_value_out)
    ELSE
       CALL WOM_Implementation(H, ISQ, K, inv_temp, params, &
            & trace_in = trace)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructSolverParameters(params)
  END SUBROUTINE WOM_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Actual implementation of WOM methods.
  SUBROUTINE WOM_Implementation(H, ISQ, K, inv_temp, params, &
       & trace_in, mu_in, energy_value_out)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: H
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: ISQ
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: K
    !> The inverse temperature for smearing in a.u.
    REAL(NTREAL), INTENT(IN) :: inv_temp
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !> The target trace.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: trace_in
    !> The target mu.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: mu_in
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !! Local Variables
    LOGICAL :: GC
    TYPE(Matrix_ps) :: ISQT, WH, IMat, RK1, RK2, K0, K1
    TYPE(Matrix_ps) :: Temp, W, A, X, KOrth
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II
    REAL(NTREAL) :: step, B_I, B_I_old, err, err2, sparsity, energy

    GC = PRESENT(mu_in)

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(IMat, H)
    CALL FillMatrixIdentity(IMat)

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL SimilarityTransform(H, ISQ, ISQT, WH, pool, &
         & threshold_in = params%threshold)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(WH, WH, &
            & params%BalancePermutation, memorypool_in = pool)
       CALL PermuteMatrix(IMat, IMat, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    !! Construct the "A" Matrix 
    IF (GC) THEN
       CALL CopyMatrix(WH, A)
       CALL IncrementMatrix(IMat, A, alpha_in = -1.0_NTREAL*mu_in)
    ELSE
       CALL CopyMatrix(WH, A)
    END IF

    !! Construct the Initial Guess
    IF (GC) THEN
       CALL CopyMatrix(IMat, W)
       CALL ScaleMatrix(W, 1.0_NTREAL/SQRT(2.0_NTREAL))
    ELSE
       CALL CopyMatrix(IMat, W)
       CALL ScaleMatrix(W, SQRT(trace_in / WH%actual_matrix_dimension))
    END IF

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF

    II = 0
    B_I = 0.0_NTREAL
    step = 1.0_NTREAL
    DO WHILE(B_I .LT. inv_temp)
       !! First order step
       step = MIN(step, inv_temp - B_I)
       CALL ComputeX(W, IMat, pool, params%threshold, X, W2_out=KOrth)
       CALL DotMatrix(WH, KOrth, energy)
       IF (GC) THEN
          CALL ComputeGCStep(X, A, pool, params%threshold, K0)
       ELSE
          CALL ComputeCStep(X, A, W, pool, params%threshold, K0)
       END IF
       II = II + 1
       CALL CopyMatrix(K0, RK1)
       CALL ScaleMatrix(RK1, step)
       CALL IncrementMatrix(W, RK1, threshold_in = params%threshold)

       !! Second Order Step
       CALL ComputeX(RK1, IMat, pool, params%threshold, X)
       IF (GC) THEN
          CALL ComputeGCStep(X, A, pool, params%threshold, K1)
       ELSE
          CALL ComputeCStep(X, A, RK1, pool, params%threshold, K1)
       END IF
       II = II + 1
       CALL CopyMatrix(W, RK2)
       CALL IncrementMatrix(K0, RK2, alpha_in = step*0.5_NTREAL, &
            & threshold_in = params%threshold)
       CALL IncrementMatrix(K1, RK2, alpha_in = step*0.5_NTREAL, &
            & threshold_in = params%threshold)

       !! Check the Error
       CALL CopyMatrix(RK1, Temp)
       CALL IncrementMatrix(RK2, Temp, alpha_in = -1.0_NTREAL, &
            & threshold_in = params%threshold)
       err = MatrixNorm(Temp)

       !! Correct Step Size as Needed
       DO WHILE (err .GT. 1.1 * params%step_thresh)
          step = step * (params%step_thresh / err)**(0.5)
          !! Update First Order
          CALL CopyMatrix(K0, RK1)
          CALL ScaleMatrix(RK1, step)
          CALL IncrementMatrix(W, RK1, threshold_in = params%threshold)
          !! Update Second Order
          CALL ComputeX(RK1, IMat, pool, params%threshold, X)
          IF (GC) THEN
             CALL ComputeGCStep(X, A, pool, params%threshold, K1)
          ELSE
             CALL ComputeCStep(X, A, RK1, pool, params%threshold, K1)
          END IF
          II = II + 1
          CALL CopyMatrix(W, RK2)
          CALL IncrementMatrix(K0, RK2, alpha_in = step*0.5_NTREAL, &
               & threshold_in = params%threshold)
          CALL IncrementMatrix(K1, RK2, alpha_in = step*0.5_NTREAL, &
               & threshold_in = params%threshold)
          !! New Error
          CALL CopyMatrix(RK1, Temp)
          CALL IncrementMatrix(RK2, Temp, alpha_in = -1.0_NTREAL, &
               & threshold_in = params%threshold)
          err = MatrixNorm(Temp)
       END DO

       !! Early Exit Criteria
       CALL CopyMatrix(RK2, Temp)
       CALL IncrementMatrix(W, Temp, alpha_in = -1.0_NTREAL, &
            & threshold_in = params%threshold)
       err2 = MatrixNorm(Temp)
       IF (err2 .LT. params%converge_diff) THEN
          CALL WriteComment("Early Exit Triggered")
          EXIT
       END IF

       !! Update
       CALL CopyMatrix(RK2, W)
       B_I_old = B_I
       B_I = B_I + step
       step = step * (params%step_thresh / err) ** (0.5)
       sparsity = REAL(GetMatrixSize(W), KIND = NTREAL) / &
            & (REAL(W%actual_matrix_dimension, KIND = NTREAL) ** 2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key = "Gradient Evaluations", VALUE = II)
          CALL EnterSubLog
          CALL WriteElement("Beta", VALUE = B_I_old)
          CALL WriteElement("Sparsity", VALUE = sparsity)
          CALL WriteElement("Energy", VALUE = energy)
          CALL WriteElement("Norm of Change", VALUE = err2)
          CALL ExitSubLog
       END IF
    END DO

    !! Form the density
    CALL MatrixMultiply(W, W, KOrth, &
         & threshold_in = params%threshold, memory_pool_in = pool)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total_Iterations", VALUE = II)
       CALL PrintMatrixInformation(W)
    END IF

    !! Compute the energy
    IF (PRESENT(energy_value_out)) THEN
       CALL DotMatrix(WH, KOrth, energy_value_out)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(KOrth, KOrth, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(KOrth, ISQT, ISQ, K, pool, &
         & threshold_in = params%threshold)

  END SUBROUTINE WOM_Implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the "X" matrix X = W [1 - W^2]
  !> Take one step for the WOM_GC algorithm. 
  SUBROUTINE ComputeX(W, I, pool, threshold, Out, W2_out)
    !> The working wave operator.
    TYPE(Matrix_ps), INTENT(IN) :: W
    !> The identity matrix.
    TYPE(Matrix_ps), INTENT(IN) :: I
    !> The memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The threshold for small values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> The result matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: Out
    !> If you want square of the wave operator
    TYPE(Matrix_ps), INTENT(INOUT), OPTIONAL :: W2_out
    !! Local matrices.
    TYPE(Matrix_ps) :: W2, Temp

    !! Build
    CALL MatrixMultiply(W, W, W2, &
         & threshold_in=threshold, memory_pool_in = pool)
    CALL CopyMatrix(W2, Temp)
    CALL ScaleMatrix(Temp, -1.0_NTREAL)
    CALL IncrementMatrix(I, Temp, threshold_in = threshold)
    CALL MatrixMultiply(W, Temp, Out, &
         & threshold_in=threshold, memory_pool_in = pool)

    IF (PRESENT(W2_out)) THEN
       CALL CopyMatrix(W2, W2_out)
    END IF

    !! Cleanup
    CALL DestructMatrix(W2)
    CALL DestructMatrix(Temp)
  END SUBROUTINE ComputeX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one step for the WOM_GC algorithm. 
  SUBROUTINE ComputeGCStep(X, A, pool, threshold, Out)
    !> The X matrix.
    TYPE(Matrix_ps), INTENT(IN) :: X
    !> H - mu*I
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The threshold for small values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> The result matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: Out

    CALL MatrixMultiply(X, A, Out, alpha_in = -0.5_NTREAL, &
         & threshold_in = threshold, memory_pool_in = pool)

  END SUBROUTINE ComputeGCStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one step for the WOM_C algorithm. 
  SUBROUTINE ComputeCStep(X, A, W, pool, threshold, Out)
    !> The X matrix.
    TYPE(Matrix_ps), INTENT(IN) :: X
    !> H
    TYPE(Matrix_ps), INTENT(IN) :: A
    !> The wave operator
    TYPE(Matrix_ps), INTENT(IN) :: W
    !> The memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The threshold for small values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> The identity matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: Out
    !! Local matrices.
    TYPE(Matrix_ps) :: XA
    REAL(NTREAL) :: num, denom

    !! Form XA
    CALL MatrixMultiply(X, A, XA, &
         & threshold_in = threshold, memory_pool_in = pool)

    !! Scaling Factor Bottom
    CALL DotMatrix(X, W, denom)
    CALL DotMatrix(W, XA, num)
    CALL CopyMatrix(X, Out)
    CALL ScaleMatrix(Out, -1.0_NTREAL * num / denom)

    !! Combine
    CALL IncrementMatrix(XA, Out)
    CALL ScaleMatrix(Out, -0.5_NTREAL)

    !! Cleanup
    CALL DestructMatrix(XA)

  END SUBROUTINE ComputeCStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorModule
