!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Wraps the load balancer module for calling from other languages.
MODULE LoadBalancerModule_wrp
  USE MatrixMemoryPoolPModule_wrp, ONLY : MatrixMemoryPool_p_wrp
  USE MatrixPSModule_wrp, ONLY : &
       & Matrix_ps_wrp
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PermuteMatrix_wrp
  PUBLIC :: UndoPermuteMatrix_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Apply a permutation to a matrix.
  SUBROUTINE PermuteMatrix_wrp(ih_mat_in, ih_mat_out, ih_permutation, &
       & ih_memorypool) bind(c, name="PermuteMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_mat_in(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat_out(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_permutation(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_memorypool(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_mat_in
    TYPE(Matrix_ps_wrp) :: h_mat_out
    TYPE(Permutation_wrp) :: h_permutation
    TYPE(MatrixMemoryPool_p_wrp) :: h_memorypool

    h_mat_in = TRANSFER(ih_mat_in,h_mat_in)
    h_mat_out = TRANSFER(ih_mat_out,h_mat_out)
    h_permutation = TRANSFER(ih_permutation,h_permutation)
    h_memorypool = TRANSFER(ih_memorypool, h_memorypool)

    CALL PermuteMatrix(h_mat_in%data, h_mat_out%data, h_permutation%data, &
         & h_memorypool%data)
  END SUBROUTINE PermuteMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Undo a permutation applied to a matrix.
  SUBROUTINE UndoPermuteMatrix_wrp(ih_mat_in, ih_mat_out, ih_permutation, &
       & ih_memorypool) bind(c, name="UndoPermuteMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_mat_in(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat_out(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_permutation(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_memorypool(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_mat_in
    TYPE(Matrix_ps_wrp) :: h_mat_out
    TYPE(Permutation_wrp) :: h_permutation
    TYPE(MatrixMemoryPool_p_wrp) :: h_memorypool

    h_mat_in = TRANSFER(ih_mat_in,h_mat_in)
    h_mat_out = TRANSFER(ih_mat_out,h_mat_out)
    h_permutation = TRANSFER(ih_permutation,h_permutation)
    h_memorypool = TRANSFER(ih_memorypool, h_memorypool)

    CALL UndoPermuteMatrix(h_mat_in%data, h_mat_out%data, h_permutation%data, &
         & h_memorypool%data)
  END SUBROUTINE UndoPermuteMatrix_wrp
END MODULE LoadBalancerModule_wrp
