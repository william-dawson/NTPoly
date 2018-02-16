!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the overlap solvers module for calling from other languages.
MODULE GeometryOptimizationModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE GeometryOptimizationModule, ONLY : PurificationExtrapolate
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PurificationExtrapolate_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !! Based on the purification algorithm in \cite niklasson2010trace .
  !! @param[in] ih_PreviousDensity to extrapolate from.
  !! @param[in] ih_Overlap the overlap matrix of the new geometry.
  !! @param[in] nel the number of electrons.
  !! @param[out] ih_NewDensity the extrapolated density.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE PurificationExtrapolate_wrp(ih_PreviousDensity, ih_Overlap, &
       & nel, ih_NewDensity, ih_solver_parameters) &
       & bind(c,name="PurificationExtrapolate_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_PreviousDensity(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_Overlap(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: nel
    INTEGER(kind=c_int), INTENT(inout) :: ih_NewDensity(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_PreviousDensity
    TYPE(DistributedSparseMatrix_wrp) :: h_Overlap
    TYPE(DistributedSparseMatrix_wrp) :: h_NewDensity
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_PreviousDensity = TRANSFER(ih_PreviousDensity,h_PreviousDensity)
    h_Overlap = TRANSFER(ih_Overlap,h_Overlap)
    h_NewDensity = TRANSFER(ih_NewDensity,h_NewDensity)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL PurificationExtrapolate(h_PreviousDensity%data, h_Overlap%data, nel, &
         & h_NewDensity%data, h_solver_parameters%data)
  END SUBROUTINE PurificationExtrapolate_wrp
END MODULE GeometryOptimizationModule_wrp
