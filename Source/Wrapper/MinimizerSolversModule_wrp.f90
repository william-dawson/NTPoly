!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the minimizer solvers module for calling from other languages.
MODULE MinimizerSolversModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule_wrp, ONLY : &
       & DistributedSparseMatrix_wrp
  USE IterativeSolversModule_wrp, ONLY : IterativeSolverParameters_wrp
  USE MinimizerSolversModule, ONLY : ConjugateGradient
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConjugateGradient_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the CG method.
  SUBROUTINE ConjugateGradient_wrp(ih_Hamiltonian, ih_InverseSquareRoot, &
       & nel, ih_Density, chemical_potential_out, ih_solver_parameters) &
       & bind(c,name="ConjugateGradient_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_Hamiltonian(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_InverseSquareRoot(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: nel
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_Density(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: chemical_potential_out
    INTEGER(kind=c_int), INTENT(IN) :: ih_solver_parameters(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_Hamiltonian
    TYPE(DistributedSparseMatrix_wrp) :: h_InverseSquareRoot
    TYPE(DistributedSparseMatrix_wrp) :: h_Density
    TYPE(IterativeSolverParameters_wrp) :: h_solver_parameters

    h_Hamiltonian = TRANSFER(ih_Hamiltonian,h_Hamiltonian)
    h_InverseSquareRoot = TRANSFER(ih_InverseSquareRoot,h_InverseSquareRoot)
    h_Density = TRANSFER(ih_Density,h_Density)
    h_solver_parameters = TRANSFER(ih_solver_parameters, h_solver_parameters)

    CALL ConjugateGradient(h_Hamiltonian%data, h_InverseSquareRoot%data, nel, &
         & h_Density%data, chemical_potential_out, h_solver_parameters%data)
  END SUBROUTINE ConjugateGradient_wrp
END MODULE MinimizerSolversModule_wrp
