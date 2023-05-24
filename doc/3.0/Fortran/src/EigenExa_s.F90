  !! Allocate Eigen Exa
  CALL InitializeEigenExa(A, nvals, PRESENT(eigenvectors_in), exa)

  !! Allocate Memory
  ALLOCATE(AD(exa%local_rows, exa%local_cols))
  AD = 0
  ALLOCATE(VD(exa%local_rows, exa%local_cols))
  VD = 0
  ALLOCATE(WD(exa%mat_dim))
  WD = 0

  !! Convert to EigenExa
#ifdef ISCOMPLEX
  CALL NTToEigen_c(A, AD, exa)
#else
  CALL NTToEigen_r(A, AD, exa)
#endif

  !! Calculate
#ifdef ISCOMPLEX
  CALL Compute_c(AD, VD, WD, exa)
#else
  CALL Compute_r(AD, VD, WD, exa)
#endif

  !! Convert Back
  IF (PRESENT(eigenvectors_in)) THEN
     CALL ConstructEmptyMatrix(eigenvectors_in, A)
#ifdef ISCOMPLEX
     CALL EigenToNT_c(VD, eigenvectors_in, params, exa)
#else
     CALL EigenToNT_r(VD, eigenvectors_in, params, exa)
#endif
  END IF

  CALL ConstructEmptyMatrix(eigenvalues, A)
  CALL ExtractEigenvalues(WD, eigenvalues, exa)

  !! Cleanup
#ifdef ISCOMPLEX
  CALL CleanUp_c(AD, VD, WD)
#else
  CALL CleanUp_r(AD, VD, WD)
#endif

  IF (params%be_verbose) THEN
     IF (PRESENT(eigenvectors_in)) THEN
        CALL PrintMatrixInformation(eigenvectors_in)
     END IF
     CALL ExitSubLog
  END IF
