!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Convert between different types of parameters.
MODULE ParameterConverterModule
  USE DataTypesModule
  USE FixedSolversModule
  USE IterativeSolversModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConvertFixedToIterative
  PUBLIC :: ConvertIterativeToFixed
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert fixed to iterative solver parameters.
  !! @param[in] fixed_parameters
  !! @param[inout] iterative_parameters
  SUBROUTINE ConvertFixedToIterative(fixed_parameters, iterative_parameters)
    TYPE(FixedSolverParameters), INTENT(IN) :: fixed_parameters
    TYPE(IterativeSolverParameters), INTENT(INOUT) :: iterative_parameters

    iterative_parameters = IterativeSolverParameters()
    iterative_parameters%threshold = fixed_parameters%threshold
    iterative_parameters%be_verbose = fixed_parameters%be_verbose
    iterative_parameters%do_load_balancing = fixed_parameters%do_load_balancing
    iterative_parameters%BalancePermutation = &
         & fixed_parameters%BalancePermutation
  END SUBROUTINE ConvertFixedToIterative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert interative to fixed solver parameters.
  !! @param[in] iterative_parameters
  !! @param[inout] fixed parameters
  SUBROUTINE ConvertIterativeToFixed(iterative_parameters, fixed_parameters)
    TYPE(IterativeSolverParameters), INTENT(IN) :: iterative_parameters
    TYPE(FixedSolverParameters), INTENT(INOUT) :: fixed_parameters

    fixed_parameters = FixedSolverParameters()
    fixed_parameters%threshold = iterative_parameters%threshold
    fixed_parameters%be_verbose = iterative_parameters%be_verbose
    fixed_parameters%do_load_balancing = iterative_parameters%do_load_balancing
    fixed_parameters%BalancePermutation = &
         & iterative_parameters%BalancePermutation
  END SUBROUTINE ConvertIterativeToFixed
END MODULE ParameterConverterModule