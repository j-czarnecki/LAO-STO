MODULE test_profiling
USE types
USE hamiltonians
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
SUBROUTINE test_profile_tba_term()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  REAL(REAL64), PARAMETER :: t_D = .5
  REAL(REAL64), PARAMETER :: t_I = .04

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_TBA_TERM")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_TBA_TERM(Hamiltonian, kx, ky, t_D, t_I, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_tba_term

SUBROUTINE test_profile_ti1_ti2()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  REAL(REAL64), PARAMETER :: eta_p = 1.
  REAL(REAL64), PARAMETER :: V_pdp = 3.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_TI1_TI2")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_TI1_TI2(Hamiltonian, kx, ky, eta_p, V_pdp, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_ti1_ti2

SUBROUTINE test_profile_h_pi()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  REAL(REAL64), PARAMETER :: eta_p = 1.
  REAL(REAL64), PARAMETER :: V_pdp = 3.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_H_PI")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_H_PI(Hamiltonian, kx, ky, eta_p, V_pdp, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_h_pi

SUBROUTINE test_profile_h_sigma()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  REAL(REAL64), PARAMETER :: eta_p = 1.
  REAL(REAL64), PARAMETER :: V_pds = 3.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_H_SIGMA")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_H_SIGMA(Hamiltonian, kx, ky, eta_p, V_pds, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_h_sigma

SUBROUTINE test_profile_rashba_hopping()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  REAL(REAL64), PARAMETER :: t_Rashba = 1.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_RASHBA_HOPPING")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky, t_Rashba, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_rashba_hopping

SUBROUTINE test_profile_gamma_sc()
  IMPLICIT NONE
  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

  !Physical variables
  COMPLEX(REAL64) :: Gamma_SC(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS) !! Superconducting energies

  REAL(REAL64), PARAMETER :: gamma_value = 1.

  INTEGER(INT32), PARAMETER :: nr_points = 2000
  INTEGER(INT32), PARAMETER :: nphi_points = 2000
  INTEGER(INT32), PARAMETER :: n_pi_3_rotations = 1
  INTEGER(INT32) :: ir, jphi
  REAL(REAL64) :: kx, ky, dr, dphi

  Gamma_SC = 0.
  !Setting fully-symmetric pairing symmetry
  Gamma_SC = gamma_value
  dr = 1./nr_points
  dphi = PI / nphi_points

  CALL timer % start("COMPUTE_SC")
  DO ir = 0, nr_points
    DO jphi = -nphi_points, nphi_points
      Hamiltonian = CMPLX(0., 0., KIND=REAL64)
      kx = ir * dr * COS(jphi * dphi)
      ky = ir * dr * SIN(jphi * dphi)
      CALL COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, sc_input % discretization)
    END DO
  END DO
  CALL timer % stop()
END SUBROUTINE test_profile_gamma_sc

END MODULE test_profiling

PROGRAM MAIN_PROFILING
USE test_profiling
IMPLICIT NONE

CALL SETUP()
CALL test_profile_tba_term()
CALL test_profile_ti1_ti2()
CALL test_profile_h_pi()
CALL test_profile_h_sigma()
CALL test_profile_rashba_hopping()
CALL test_profile_gamma_sc()

END PROGRAM MAIN_PROFILING
