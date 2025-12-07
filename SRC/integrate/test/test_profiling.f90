MODULE test_profiling
USE types
USE local_integrand
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
SUBROUTINE test_profile_accumulate_nearest_neighbors_delta()
  IMPLICIT NONE
  COMPLEX(REAL64):: Delta(sc_input % discretization % ORBITALS, &
                         & N_ALL_NEIGHBOURS, &
                         & SPINS, &
                         & SPINS, &
                         & sc_input % discretization % derived % LAYER_COUPLINGS) !! Accumulator for integrand
  REAL(REAL64) :: J_tensor(SPINS, SPINS, SPINS, SPINS)
  COMPLEX(REAL64) :: U(sc_input % discretization % derived % DIM, &
                      & sc_input % discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64) :: Energies(sc_input % discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), PARAMETER :: T = 0.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  Delta = 0.
  J_tensor = 1.
  U = 1.
  Energies = -1.

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("ACCUMULATE_NEAREST_NEIGHBORS_DELTA")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL ACCUMULATE_NEAREST_NEIGHBOURS_DELTA(Delta, J_tensor, U, Energies, kx, ky, sc_input % discretization, T)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_accumulate_nearest_neighbors_delta

SUBROUTINE test_profile_accumulate_next_neighbors_delta()
  IMPLICIT NONE
  COMPLEX(REAL64):: Delta(sc_input % discretization % ORBITALS, &
                         & N_ALL_NEIGHBOURS, &
                         & SPINS, &
                         & SPINS, &
                         & sc_input % discretization % derived % LAYER_COUPLINGS) !! Accumulator for integrand
  REAL(REAL64) :: J_tensor(SPINS, SPINS, SPINS, SPINS)
  COMPLEX(REAL64) :: U(sc_input % discretization % derived % DIM, &
                      & sc_input % discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64) :: Energies(sc_input % discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), PARAMETER :: T = 0.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  Delta = 0.
  J_tensor = 1.
  U = 1.
  Energies = -1.

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("ACCUMULATE_NEXT_NEIGHBORS_DELTA")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL ACCUMULATE_NEXT_NEIGHBOURS_DELTA(Delta, J_tensor, U, Energies, kx, ky, sc_input % discretization, T)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_accumulate_next_neighbors_delta

END MODULE test_profiling

PROGRAM MAIN_PROFILING
USE test_profiling
IMPLICIT NONE

CALL SETUP()
CALL test_profile_accumulate_nearest_neighbors_delta()
CALL test_profile_accumulate_next_neighbors_delta()
END PROGRAM MAIN_PROFILING
