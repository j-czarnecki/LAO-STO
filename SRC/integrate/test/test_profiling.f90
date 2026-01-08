MODULE test_profiling
USE types
USE local_integrand
USE utilities
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
IMPLICIT NONE

!------------------------------------------------------------------------------
!-------------------------------- TYPES ---------------------------------------
!------------------------------------------------------------------------------
TYPE :: timer_t
  REAL(REAL64) :: t_start
  REAL(REAL64) :: t_end
  REAL(REAL64) :: t_elapsed
CONTAINS
  PROCEDURE :: start
  PROCEDURE :: stop
END TYPE timer_t

!------------------------------------------------------------------------------
!----------------------------- LOCAL VARIAVBLES -------------------------------
!------------------------------------------------------------------------------
REAL(REAL64), PARAMETER :: s = 0.5
TYPE(sc_input_params_t) :: sc_input
TYPE(timer_t) :: timer

CONTAINS
!------------------------------------------------------------------------------
!---------------------------- SETUP/TEARDOWN ----------------------------------
!------------------------------------------------------------------------------
SUBROUTINE setup()
  INTEGER(INT32) :: sublats = 2
  INTEGER(INT32) :: n_subbands = 2
  CALL SET_HAMILTONIAN_PARAMS(sublats, n_subbands, sc_input % discretization)
END SUBROUTINE setup

SUBROUTINE start(this, function_name)
  IMPLICIT NONE
  CLASS(timer_t), INTENT(INOUT) :: this
  CHARACTER(LEN=*), INTENT(IN) :: function_name
  INTEGER(INT32) :: count, rate
  WRITE (*, *) "Measuring function: ", function_name
  CALL SYSTEM_CLOCK(count, rate)
  this % t_start = REAL(count, REAL64) / REAL(rate, REAL64)
END SUBROUTINE start

SUBROUTINE stop(this)
  IMPLICIT NONE
  CLASS(timer_t), INTENT(INOUT) :: this
  INTEGER(INT32) :: count, rate
  CALL SYSTEM_CLOCK(count, rate)
  this % t_end = REAL(count, REAL64) / REAL(rate, REAL64)
  this % t_elapsed = this % t_end - this % t_start
  WRITE (*, *) "Time elapsed: ", this % t_elapsed
END SUBROUTINE stop

!---------------------------------------------------------------------
!----------------------------- TESTS ---------------------------------
!---------------------------------------------------------------------
SUBROUTINE test_profile_accumulate_delta_real_space()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Delta(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K) !! Accumulator for integrand
  REAL(REAL64) :: J_tensor(sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K) !! Interaction tensor
  REAL(REAL64) :: J_tensor_matrix(sc_input % discretization % derived % DIM_POSITIVE_K**2, &
                                & sc_input % discretization % derived % DIM_POSITIVE_K**2) !! Interaction tensor
  REAL(REAL64), ALLOCATABLE :: Values_j_tensor(:) !! Interaction tensor
  INTEGER(INT32), ALLOCATABLE :: Column_indices_j_tensor(:) !! Interaction tensor
  INTEGER(INT32), ALLOCATABLE :: Row_indices_j_tensor(:) !! Interaction tensor
  COMPLEX(REAL64) :: U(sc_input % discretization % derived % DIM, &
                      & sc_input % discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64) :: Energies(sc_input % discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), PARAMETER :: T = 0.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi, nonzero_elems, i
  REAL(REAL64) :: kx, ky, dr, dphi

  Delta = 0.
  J_tensor = 0.
  DO i = 1, sc_input % discretization % derived % DIM_POSITIVE_K
    J_tensor(i, i, i, i) = 1.
  END DO
  U = 1.
  Energies = -1.

  CALL MATRICIZE_INTERACTION_TENSOR(J_tensor, sc_input % discretization % derived % DIM_POSITIVE_K, J_tensor_matrix)
  nonzero_elems = calculate_number_of_nonzero_elements(J_tensor_matrix)
  WRITE (*, *) "Nonzero elements: ", nonzero_elems

  ALLOCATE (Values_j_tensor(nonzero_elems))
  ALLOCATE (Column_indices_j_tensor(nonzero_elems))
  ALLOCATE (Row_indices_j_tensor(sc_input % discretization % derived % DIM_POSITIVE_K**2 + 1))
  CALL SAVE_SPARSE_MATRIX_IN_CRS(J_tensor_matrix, Values_j_tensor, Column_indices_j_tensor, Row_indices_j_tensor, &
    & nonzero_elems, sc_input % discretization % derived % DIM_POSITIVE_K**2)

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("ACCUMULATE_NEAREST_NEIGHBORS_DELTA")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL ACCUMULATE_DELTA_REAL_SPACE(Delta, Values_j_tensor, Column_indices_j_tensor, Row_indices_j_tensor, &
        & U, Energies, kx, ky, sc_input % discretization, nonzero_elems, T)
    END DO
  END DO
  CALL timer % stop()

  DEALLOCATE (Values_j_tensor)
  DEALLOCATE (Column_indices_j_tensor)
  DEALLOCATE (Row_indices_j_tensor)
END SUBROUTINE test_profile_accumulate_delta_real_space

SUBROUTINE test_profile_accumulate_delta_k_space()
  COMPLEX(REAL64) :: Delta(sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K) !! Accumulator for integrand
  REAL(REAL64) :: J_tensor(sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K, &
                         & sc_input % discretization % derived % DIM_POSITIVE_K) !! Interaction tensor
  REAL(REAL64) :: J_tensor_matrix(sc_input % discretization % derived % DIM_POSITIVE_K**2, &
                                & sc_input % discretization % derived % DIM_POSITIVE_K**2) !! Interaction tensor
  REAL(REAL64), ALLOCATABLE :: Values_j_tensor(:) !! Interaction tensor
  INTEGER(INT32), ALLOCATABLE :: Column_indices_j_tensor(:) !! Interaction tensor
  INTEGER(INT32), ALLOCATABLE :: Row_indices_j_tensor(:) !! Interaction tensor

  COMPLEX(REAL64) :: U(sc_input % discretization % derived % DIM, &
                     & sc_input % discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64) :: Energies(sc_input % discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64) :: kx, ky !! Wavevector coordinates
  REAL(REAL64) :: T !! Temperature

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: dr, dphi

  INTEGER(INT32) :: i, j, nonzero_elems

  Delta = 0.
  DO i = 1, sc_input % discretization % derived % DIM_POSITIVE_K
    j = MOD(i, sc_input % discretization % derived % DIM_POSITIVE_K) + 1
    J_tensor(i, j, i, j) = 1.
    J_tensor(j, i, j, i) = 1.
    J_tensor(i, i, j, j) = 1.
    J_tensor(i, i, i, i) = 1.
  END DO
  U = 1.
  Energies = -1.

  !! Prepare for function
  CALL MATRICIZE_INTERACTION_TENSOR(J_tensor, sc_input % discretization % derived % DIM_POSITIVE_K, J_tensor_matrix)

  nonzero_elems = calculate_number_of_nonzero_elements(J_tensor_matrix)
  WRITE (*, *) "Nonzero elements: ", nonzero_elems

  ALLOCATE (Values_j_tensor(nonzero_elems))
  ALLOCATE (Column_indices_j_tensor(nonzero_elems))
  ALLOCATE (Row_indices_j_tensor(sc_input % discretization % derived % DIM_POSITIVE_K**2 + 1))

  CALL SAVE_SPARSE_MATRIX_IN_CRS(J_tensor_matrix, Values_j_tensor, Column_indices_j_tensor, Row_indices_j_tensor, &
    & nonzero_elems, sc_input % discretization % derived % DIM_POSITIVE_K**2)

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("ACCUMULATE_DELTA_K_SPACE")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL ACCUMULATE_DELTA_K_SPACE(Delta, Values_j_tensor, Column_indices_j_tensor, Row_indices_j_tensor, &
        & U, Energies, kx, ky, sc_input % discretization, nonzero_elems, T)
    END DO
  END DO
  CALL timer % stop()

  DEALLOCATE (Values_j_tensor)
  DEALLOCATE (Column_indices_j_tensor)
  DEALLOCATE (Row_indices_j_tensor)

END SUBROUTINE test_profile_accumulate_delta_k_space

END MODULE test_profiling

PROGRAM MAIN_PROFILING
USE test_profiling
IMPLICIT NONE

CALL SETUP()
CALL test_profile_accumulate_delta_real_space()
CALL test_profile_accumulate_delta_k_space()
END PROGRAM MAIN_PROFILING
