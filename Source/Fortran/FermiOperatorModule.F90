!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Density Matrix Using the Fermi Operator Expansion
MODULE FermiOperatorModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE EigenSolversModule, ONLY : EigenDecomposition
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteHeader, &
       & EnterSubLog, ExitSubLog, WriteListElement
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, SimilarityTransform, &
       & IncrementMatrix, MatrixNorm, ScaleMatrix, DotMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & FillMatrixFromTripletList, GetMatrixTripletList, &
       & TransposeMatrix, ConjugateMatrix, DestructMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, CopyMatrix, GetMatrixSize
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & PrintParameters, DestructSolverParameters, CopySolverParameters, &
       & ConstructSolverParameters
  USE TripletListModule, ONLY : TripletList_r, DestructTripletList
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeDenseFOE
  PUBLIC :: WOM_GC
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
    REAL(NTREAL) :: chemical_potential, energy_value
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: eigs, occ
    REAL(NTREAL) :: sval, sv, occ_temp
    REAL(NTREAL) :: left, right, homo, lumo
    INTEGER :: num_eigs
    INTEGER :: II, JJ
    INTEGER :: ierr

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
          CALL WriteElement(key="Method", VALUE="Dense FOE")
          CALL WriteElement(key="Inverse Temperature", VALUE=inv_temp)
       ELSE
          CALL WriteElement(key="Method", VALUE="Dense Step Function")
       END IF
       CALL PrintParameters(params)
    END IF

    !! Compute the working hamiltonian.
    CALL TransposeMatrix(ISQ, ISQT)
    CALL MatrixMultiply(ISQ, H, Temp, &
         & threshold_in=params%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(Temp, ISQT, WH, &
         & threshold_in=params%threshold, memory_pool_in=pool)

    !! Perform the eigendecomposition
    CALL EigenDecomposition(WH, vals, &
         & eigenvectors_in=vecs, solver_parameters_in=params)

    !! Gather the eigenvalues on to every process
    CALL GetMatrixTripletList(vals, tlist)
    num_eigs = H%actual_matrix_dimension
    ALLOCATE(eigs(num_eigs))
    eigs = 0
    DO II = 1, tlist%CurrentSize
       eigs(tlist%DATA(II)%index_column) = tlist%DATA(II)%point_value
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, eigs, num_eigs, MPINTREAL, &
         & MPI_SUM, H%process_grid%within_slice_comm, ierr)

    !! Compute MU By Bisection
    IF (do_smearing) THEN
       ALLOCATE(occ(num_eigs))
       left = MINVAL(eigs)
       right = MAXVAL(eigs)
       DO JJ = 1, 10*params%max_iterations 
          chemical_potential = left + (right - left)/2 
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
       CALL WriteElement(key="Potential", VALUE=chemical_potential)
       CALL WriteElement(key="Iterations", VALUE=JJ)
       CALL ExitSubLog
    END IF

    !! Map
    energy_value = 0.0_NTREAL
    DO II = 1, tlist%CurrentSize
       IF (.NOT. do_smearing) THEN
          IF (tlist%DATA(II)%index_column .LE. FLOOR(trace)) THEN
             energy_value = energy_value + tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = 1.0_NTREAL
          ELSE IF (tlist%DATA(II)%index_column .EQ. CEILING(trace)) THEN
             occ_temp = CEILING(trace) - trace
             energy_value = energy_value + &
                  & occ_temp * tlist%DATA(II)%point_value
             tlist%DATA(II)%point_value = occ_temp
          ELSE
             tlist%DATA(II)%point_value = 0.0_NTREAL
          ENDIF
       ELSE
          sval = tlist%DATA(II)%point_value - chemical_potential
          ! occ_temp = 0.5_NTREAL * (1.0_NTREAL - ERF(inv_temp * sval))
          occ_temp = 1.0_NTREAL / (1.0_NTREAL + EXP(inv_temp * sval))
          energy_value = energy_value + occ_temp * tlist%DATA(II)%point_value
          tlist%DATA(II)%point_value = occ_temp
       END IF
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, energy_value, 1, MPINTREAL, MPI_SUM, &
         & H%process_grid%within_slice_comm, ierr)

    !! Fill
    CALL ConstructEmptyMatrix(vals, H)
    CALL FillMatrixFromTripletList(vals, tlist, preduplicated_in=.TRUE.)

    !! Multiply Back Together
    CALL MatrixMultiply(vecs, vals, temp, threshold_in=params%threshold)
    CALL TransposeMatrix(vecs, vecsT)
    CALL ConjugateMatrix(vecsT)
    CALL MatrixMultiply(temp, vecsT, WD, &
         & threshold_in=params%threshold)

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(ISQT, WD, Temp, &
         & threshold_in=params%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(Temp, ISQ, K, &
         & threshold_in=params%threshold, memory_pool_in=pool)

    !! Optional out variables.
    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = 2.0_NTREAL * energy_value
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
    !! Local Variables
    TYPE(Matrix_ps) :: ISQT, WH, IMat, RK1, RK2, K0, K1
    TYPE(Matrix_ps) :: Temp, W, X, KOrth
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: II
    REAL(NTREAL) :: step, B_I, err, sparsity

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
       CALL WriteElement(key="Method", VALUE="WOM_GC")
       CALL WriteElement(key="Inverse Temperature", VALUE=inv_temp)
       CALL WriteElement(key="Chemical Potential", VALUE=chemical_potential)
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
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

    !! Construct the "X" Matrix 
    CALL CopyMatrix(WH, X)
    CALL IncrementMatrix(IMat, X, alpha_in=-1.0_NTREAL*chemical_potential)

    !! Construct the Initial Guess
    CALL CopyMatrix(IMat, W)
    CALL ScaleMatrix(W, 1.0_NTREAL/SQRT(2.0_NTREAL))

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
       step = min(step, inv_temp - B_I)
       CALL ComputeGCStep(W, IMat, X, pool, params%threshold, K0)
       II = II + 1
       CALL CopyMatrix(K0, RK1)
       CALL ScaleMatrix(RK1, step)
       CALL IncrementMatrix(W, RK1, threshold_in=params%threshold)

       !! Second Order Step
       CALL ComputeGCStep(RK1, IMat, X, pool, params%threshold, K1)
       II = II + 1
       CALL CopyMatrix(W, RK2)
       CALL IncrementMatrix(K0, RK2, alpha_in=step*0.5_NTREAL, &
            & threshold_in=params%threshold)
       CALL IncrementMatrix(K1, RK2, alpha_in=step*0.5_NTREAL, &
            & threshold_in=params%threshold)

       !! Check the Error
       CALL CopyMatrix(RK1, Temp)
       CALL IncrementMatrix(RK2, Temp, alpha_in=-1.0_NTREAL, &
            & threshold_in=params%threshold)
       err = MatrixNorm(Temp)

       !! Correct Step Size as Needed
       DO WHILE (err .GT. 1.1 * params%step_thresh)
          step = step * (params%step_thresh / err)**(0.5)
          !! Update First Order
          CALL CopyMatrix(K0, RK1)
          CALL ScaleMatrix(RK1, step)
          CALL IncrementMatrix(W, RK1, threshold_in=params%threshold)
          !! Update Second Order
          CALL ComputeGCStep(RK1, IMat, X, pool, params%threshold, K1)
          II = II + 1
          CALL CopyMatrix(W, RK2)
          CALL IncrementMatrix(K0, RK2, alpha_in=step*0.5_NTREAL, &
               & threshold_in=params%threshold)
          CALL IncrementMatrix(K1, RK2, alpha_in=step*0.5_NTREAL, &
               & threshold_in=params%threshold)
          !! New Error
          CALL CopyMatrix(RK1, Temp)
          CALL IncrementMatrix(RK2, Temp, alpha_in=-1.0_NTREAL, &
               & threshold_in=params%threshold)
          err = MatrixNorm(Temp)
       END DO

       !! Update
       CALL CopyMatrix(RK2, W)
       B_I = B_I + step
       step = step * (params%step_thresh / err)**(0.5)
       sparsity = REAL(GetMatrixSize(W),KIND=NTREAL) / &
            & (REAL(W%actual_matrix_dimension,KIND=NTREAL)**2)

       IF (params%be_verbose) THEN
          CALL WriteListElement(key="GC Steps", VALUE=II)
          CALL EnterSubLog
          CALL WriteElement("Beta", VALUE=B_I)
          CALL WriteElement("Sparsity", VALUE=sparsity)
          CALL ExitSubLog
       END IF
    END DO

    !! Form the density
    CALL MatrixMultiply(W, W, KOrth, &
         & threshold_in=params%threshold, memory_pool_in=pool)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=II)
       CALL PrintMatrixInformation(W)
    END IF

    !! Compute the energy
    CALL DotMatrix(WH, KOrth, energy_value_out)

     !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(KOrth, KOrth, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL SimilarityTransform(KOrth, ISQT, ISQ, K, pool, &
         & threshold_in=params%threshold)

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
    CALL DestructMatrix(ISQT)
    CALL DestructMatrix(WH)
    CALL DestructMatrix(IMat)
    CALL DestructMatrix(RK1)
    CALL DestructMatrix(RK2)
    CALL DestructMatrix(K0)
    CALL DestructMatrix(K1)
    CALL DestructMatrix(Temp)
    CALL DestructMatrix(W)
    CALL DestructMatrix(X)
    CALL DestructMatrix(KOrth)
    CALL DestructSolverParameters(params)
  END SUBROUTINE WOM_GC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take one step for the WOM_GC algorithm. 
  SUBROUTINE ComputeGCStep(W, I, X, pool, threshold, Out)
    !> The working wave operator.
    TYPE(Matrix_ps), INTENT(IN) :: W
    !> The identity matrix.
    TYPE(Matrix_ps), INTENT(IN) :: I
    !> H - mu*I
    TYPE(Matrix_ps), INTENT(IN) :: X
    !> The memory pool.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: pool
    !> The threshold for small values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> The identity matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: Out
    !! Local matrices.
    TYPE(Matrix_ps) :: W2, Temp, Temp2

    !! Build
    CALL MatrixMultiply(W, W, W2, &
         & threshold_in=threshold, memory_pool_in=pool)
    CALL CopyMatrix(W2, Temp)
    CALL ScaleMatrix(Temp, -1.0_NTREAL)
    CALL IncrementMatrix(I, Temp, threshold_in=threshold)
    CALL MatrixMultiply(W, Temp, Temp2, &
         & threshold_in=threshold, memory_pool_in=pool)
    CALL MatrixMultiply(Temp2, X, Out, alpha_in=-0.5_NTREAL, &
         & threshold_in=threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(W2)
    CALL DestructMatrix(Temp)
    CALL DestructMatrix(Temp2)
  END SUBROUTINE ComputeGCStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FermiOperatorModule
