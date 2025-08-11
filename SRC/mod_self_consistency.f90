!! This file is part of LAO-STO.
!!
!! Copyright (C) 2025 Julian Czarnecki
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! If you use this code for scientific research, please cite:
!! J. Czarnecki et. al.,
!! "Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface",
!! arXiv:2508.05075 (2025).
!! https://arxiv.org/abs/2508.05075

#include "macros_def.f90"

MODULE mod_self_consistency
USE mod_parameters
USE mod_reader
USE mod_logger

IMPLICIT NONE
CONTAINS

SUBROUTINE GET_GAMMAS_FROM_DELTAS(Gamma_SC, Delta)
  COMPLEX*16, INTENT(INOUT) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  COMPLEX*16, INTENT(IN) :: Delta(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)

  INTEGER*4 :: band, spin, n, lat, orb, orb_prime, band_prime

  !Gamma calculation
  DO band = 1, SUBBANDS
    DO spin = 1, 2   !Loop over spin coupling up-down or down-up
      DO n = 1, N_NEIGHBOURS !Loop over neighbours
        DO lat = 1, LAYER_COUPLINGS !Loop over LAYER_COUPLINGS to include Ti1-Ti2 and Ti2-Ti1 coupling
          DO orb = 1, ORBITALS

            Gamma_SC(orb, n, spin, lat, band) = -0.5 * J_SC * Delta(orb, n, spin, lat, band)
            ! To get same critical temperature in all orbitals
            DO orb_prime = 1, ORBITALS
              IF (orb .NE. orb_prime) THEN
                Gamma_SC(orb, n, spin, lat, band) = Gamma_SC(orb, n, spin, lat, band) - 0.5 * J_SC_PRIME * Delta(orb_prime, n, spin, lat, band)
              END IF
            END DO
            !To get same critical temperature in all subbands
            DO band_prime = 1, SUBBANDS
              IF (band .NE. band_prime) THEN
                DO orb_prime = 1, ORBITALS
                  Gamma_SC(orb, n, spin, lat, band) = Gamma_SC(orb, n, spin, lat, band) - 0.5 * J_SC_PRIME / 10.0d0 * Delta(orb_prime, n, spin, lat, band_prime)
                END DO
              END IF
            END DO

          END DO
        END DO
      END DO
    END DO
  END DO
  !Gamma for next nearest neighbours pairing
  DO band = 1, SUBBANDS
    DO spin = 1, 2   !Loop over spin coupling up-down or down-up
      DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS !Loop over next nearest neighbours
        DO lat = 1, SUBLATTICES !Loop over sublattices coupling, since for next-to-nearest neighbours it has intra-layer character
          DO orb = 1, ORBITALS

            Gamma_SC(orb, n, spin, lat, band) = -0.5 * J_SC_NNN * Delta(orb, n, spin, lat, band)
            ! To get same critical temperature in all orbitals
            DO orb_prime = 1, ORBITALS
              IF (orb .NE. orb_prime) THEN
                ! CHECK WHETHER THIS IS 0.25 OR 0.5
                Gamma_SC(orb, n, spin, lat, band) = Gamma_SC(orb, n, spin, lat, band) - 0.5 * J_SC_PRIME_NNN * Delta(orb_prime, n, spin, lat, band)
              END IF
            END DO
            !To get same critical temperature in all subbands
            DO band_prime = 1, SUBBANDS
              IF (band .NE. band_prime) THEN
                DO orb_prime = 1, ORBITALS
                  Gamma_SC(orb, n, spin, lat, band) = Gamma_SC(orb, n, spin, lat, band) - 0.5 * J_SC_PRIME_NNN / 10.0d0 * Delta(orb_prime, n, spin, lat, band_prime)
                END DO
              END IF
            END DO

          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE GET_GAMMAS_FROM_DELTAS

SUBROUTINE CHECK_CONVERGENCE(sc_flag, Gamma_old, Gamma_new, Charge_dens, Charge_dens_new, gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error, sc_iter)
  LOGICAL, INTENT(INOUT) :: sc_flag
  COMPLEX*16, INTENT(IN) :: Gamma_old(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  COMPLEX*16, INTENT(IN) :: Gamma_new(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS)
  REAL*8, INTENT(IN) :: Charge_dens_new(DIM_POSITIVE_K, SUBBANDS)
  REAL*8, INTENT(INOUT) :: gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error
  INTEGER*4, INTENT(IN) :: sc_iter
  INTEGER*4 :: band, spin, orb, n, lat
  REAL*8 :: gamma_error, charge_error

  gamma_max_error_prev = gamma_max_error
  charge_max_error_prev = charge_max_error
  gamma_max_error = 0.
  charge_max_error = 0.
  !Here we check whether convergence was reached
  sc_flag = .TRUE.
  DO band = 1, SUBBANDS
    DO spin = 1, 2
      DO orb = 1, ORBITALS
        DO n = 1, N_NEIGHBOURS
          DO lat = 1, LAYER_COUPLINGS
            !It should be considered whether relative or absolute error must be checked
            gamma_error = ABS(ABS(Gamma_new(orb, n, spin, lat, band)) - ABS(Gamma_old(orb, n, spin, lat, band)))
            !Gamma convergence checking
            IF (gamma_error > gamma_eps_convergence) THEN
              sc_flag = .FALSE.
              !EXIT !Maybe go to???
            END IF

            !Find biggest error in current iteration
            IF (gamma_error > gamma_max_error) gamma_max_error = gamma_error

          END DO
        END DO
        DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
          DO lat = 1, SUBLATTICES
            !It should be considered whether relative or absolute error must be checked
            gamma_error = ABS(ABS(Gamma_new(orb, n, spin, lat, band)) - ABS(Gamma_old(orb, n, spin, lat, band)))
            !Gamma convergence checking
            IF (gamma_error > gamma_eps_convergence) THEN
              sc_flag = .FALSE.
              !EXIT !Maybe go to???
            END IF

            !Find biggest error in current iteration
            IF (gamma_error > gamma_max_error) gamma_max_error = gamma_error

          END DO
        END DO
      END DO
    END DO
  END DO
  !Change Broyden mixing parameter if simulation diverges between iterations
  !To avoid oscillations near convergence
  IF (ABS(gamma_max_error) > ABS(gamma_max_error_prev) .AND. (sc_iter > 1)) THEN
    !PRINT*, "Adapted sc_alpha = ", sc_alpha
    WRITE (log_string, '(a, E15.5)') "Divergent iteration. Adapted sc_alpha: ", sc_alpha
    LOG_ABNORMAL(log_string)
    sc_alpha = sc_alpha * sc_alpha_adapt
  END IF

  DO band = 1, SUBBANDS
    DO n = 1, DIM_POSITIVE_K
      charge_error = ABS(Charge_dens(n, band) - Charge_dens_new(n, band))
      IF (charge_error > charge_eps_convergence) THEN
        sc_flag = .FALSE.
      END IF

      IF (charge_error > charge_max_error) charge_max_error = charge_error
    END DO
  END DO

END SUBROUTINE CHECK_CONVERGENCE

SUBROUTINE FLATTEN_FOR_BROYDEN(Gamma_old, Gamma_new, Charge_dens, Charge_dens_new, Broyden_vector_new, Broyden_vector, broyden_length)
  COMPLEX*16, INTENT(IN) :: Gamma_old(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  COMPLEX*16, INTENT(IN) :: Gamma_new(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS)
  REAL*8, INTENT(IN) :: Charge_dens_new(DIM_POSITIVE_K, SUBBANDS)
  INTEGER*4, INTENT(IN) :: broyden_length
  REAL*8, INTENT(OUT) :: Broyden_vector_new(broyden_length)
  REAL*8, INTENT(OUT) :: Broyden_vector(broyden_length)

  INTEGER*4 :: band, spin, orb, n, lat, broyden_index

  broyden_index = 1
  DO band = 1, SUBBANDS
    DO spin = 1, 2
      DO orb = 1, ORBITALS
        DO n = 1, N_NEIGHBOURS
          DO lat = 1, LAYER_COUPLINGS
            Broyden_vector(broyden_index) = REAL(Gamma_old(orb, n, spin, lat, band))
            Broyden_vector_new(broyden_index) = REAL(Gamma_new(orb, n, spin, lat, band))
            broyden_index = broyden_index + 1

            Broyden_vector(broyden_index) = AIMAG(Gamma_old(orb, n, spin, lat, band))
            Broyden_vector_new(broyden_index) = AIMAG(Gamma_new(orb, n, spin, lat, band))
            broyden_index = broyden_index + 1
          END DO
        END DO
        DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
          DO lat = 1, SUBLATTICES
            Broyden_vector(broyden_index) = REAL(Gamma_old(orb, n, spin, lat, band))
            Broyden_vector_new(broyden_index) = REAL(Gamma_new(orb, n, spin, lat, band))
            broyden_index = broyden_index + 1

            Broyden_vector(broyden_index) = AIMAG(Gamma_old(orb, n, spin, lat, band))
            Broyden_vector_new(broyden_index) = AIMAG(Gamma_new(orb, n, spin, lat, band))
            broyden_index = broyden_index + 1
          END DO
        END DO

      END DO
    END DO
  END DO
  !Must be +1!!!
  DO band = 1, SUBBANDS
    DO n = 1, DIM_POSITIVE_K
      Broyden_vector(broyden_index) = Charge_dens(n, band)
      Broyden_vector_new(broyden_index) = Charge_dens_new(n, band)
      broyden_index = broyden_index + 1
    END DO
  END DO
  !Sanity check
  IF (broyden_index - 1 /= broyden_length) THEN
    WRITE (log_string, *) 'Broyden index - 1 /= delta_real_elems', broyden_index - 1, broyden_length
    LOG_ERROR(log_string)
  END IF
END SUBROUTINE FLATTEN_FOR_BROYDEN

SUBROUTINE RESHAPE_FROM_BROYDEN(Gamma_SC, Charge_dens, Broyden_vector, broyden_length)
  COMPLEX*16, INTENT(OUT) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
  REAL*8, INTENT(OUT) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS)
  INTEGER*4, INTENT(IN) :: broyden_length
  REAL*8, INTENT(IN) :: Broyden_vector(broyden_length)

  INTEGER*4 :: band, spin, orb, n, lat, broyden_index
  broyden_index = 1
  DO band = 1, SUBBANDS
    DO spin = 1, 2
      DO orb = 1, ORBITALS
        DO n = 1, N_NEIGHBOURS
          DO lat = 1, LAYER_COUPLINGS
            Gamma_SC(orb, n, spin, lat, band) = DCMPLX(Broyden_vector(broyden_index), Broyden_vector(broyden_index + 1))
            broyden_index = broyden_index + 2
          END DO
        END DO
        DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
          DO lat = 1, SUBLATTICES
            Gamma_SC(orb, n, spin, lat, band) = DCMPLX(Broyden_vector(broyden_index), Broyden_vector(broyden_index + 1))
            broyden_index = broyden_index + 2
          END DO
        END DO
      END DO
    END DO
  END DO
  DO band = 1, SUBBANDS
    DO n = 1, DIM_POSITIVE_K
      Charge_dens(n, band) = Broyden_vector(broyden_index)
      broyden_index = broyden_index + 1
    END DO
  END DO
  !Sanity check
  IF (broyden_index - 1 /= broyden_length) THEN
    WRITE (log_string, *) 'Broyden index - 1 /= delta_real_elems', broyden_index - 1, broyden_length
    LOG_ERROR(log_string)
  END IF

END SUBROUTINE

END MODULE mod_self_consistency
