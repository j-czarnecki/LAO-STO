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

MODULE self_consistency
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
USE reader
USE logger
USE types
USE writers
IMPLICIT NONE
CONTAINS

SUBROUTINE SET_GAMMA_INITIAL(Gamma_SC, J_tensor_crs, gamma_start_nearest, gamma_start_next, discretization)
  !! This subroutine sets initial values of superconducting couplings
  !! based on the energy tensors i.e. sets values relevant to the non-zero elements of tensors
  !! to gamma_start_nearest(next) and others to zero.
  !! TODO: This should eventually enable all possible pairing and randomize the initial condition.
  !! This way, it could be run several times to find the lowest-energy solution.
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization !! Discretization parameters used in calculations
  TYPE(crs_matrix_t), INTENT(IN) :: J_tensor_crs !! Interaction tensor for superconducting coupling
  REAL(REAL64), INTENT(IN) :: gamma_start_nearest, gamma_start_next !! Starting value of superconducting coupling for nearest neighbours and next neighbours
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(OUT) :: Gamma_SC(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#else
  COMPLEX(REAL64), INTENT(OUT) :: Gamma_SC(discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#endif
  INTEGER(INT32) :: spin1, spin2, spin3, spin4, spin_opposite
  LOGICAL :: set_to_nonzero_nearest_opposite_spin, set_to_nonzero_next_opposite_spin
  LOGICAL :: set_to_nonzero_nearest_same_spin, set_to_nonzero_next_same_spin
  REAL(REAL64), PARAMETER :: eps = 1e-9 !! To compare reals
  REAL(REAL64), PARAMETER :: spin_offset_fraction = 0.05
  REAL(REAL64) :: spin_sign
  INTEGER(INT32) :: i
  INTEGER(INT32) :: Dematricized_indices(2)
  INTEGER(INT32) :: last_crs_row_index
  Gamma_SC = CMPLX(0.0, 0.0, KIND=REAL64)

  DO i = 1, discretization % derived % DIM_POSITIVE_K**2
    !! If he row pointer are the same, it means that i-th row consists of zeros only
    IF (J_tensor_crs % Row_indices(i + 1) .NE. J_tensor_crs % Row_indices(i)) THEN
      Dematricized_indices = get_dematricized_indeces(i, discretization % derived % DIM_POSITIVE_K)
#ifndef BAND_BASIS
      Gamma_SC(:N_NEAREST_NEIGHBOURS, Dematricized_indices(1), Dematricized_indices(2), :) = gamma_start_nearest
      Gamma_SC(N_NEAREST_NEIGHBOURS + 1:N_ALL_NEIGHBOURS, Dematricized_indices(1), Dematricized_indices(2), :) = gamma_start_next
#else
      Gamma_SC(Dematricized_indices(1), Dematricized_indices(2), :) = gamma_start_nearest * (-1.0)**(Dematricized_indices(1))
#endif
    END IF
  END DO

  !Printing initial values
  CALL PRINT_GAMMA(Gamma_SC, "Gamma_SC_initial", discretization)

END SUBROUTINE SET_GAMMA_INITIAL

SUBROUTINE GET_GAMMAS_FROM_DELTAS(Gamma_SC, Delta, discretization, nearest_interorb_multiplier, next_interorb_multiplier)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  COMPLEX(REAL64), INTENT(INOUT) :: Gamma_SC(discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, discretization % derived % LAYER_COUPLINGS, discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(IN) :: Delta(discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, discretization % derived % LAYER_COUPLINGS, discretization % SUBBANDS)
  REAL(REAL64), INTENT(IN) :: nearest_interorb_multiplier, next_interorb_multiplier

  INTEGER(INT32) :: band, spin1, spin2, n, lat, orb, orb_prime, band_prime

  !Gamma calculation
  DO band = 1, discretization % SUBBANDS
    DO spin1 = 1, SPINS   !Loop over spin coupling up-down or down-up
      DO spin2 = 1, SPINS
        DO n = 1, N_NEIGHBOURS !Loop over neighbours
          DO lat = 1, discretization % derived % LAYER_COUPLINGS !Loop over LAYER_COUPLINGS to include Ti1-Ti2 and Ti2-Ti1 coupling
            DO orb = 1, discretization % ORBITALS

              Gamma_SC(orb, n, spin1, spin2, lat, band) = -0.5 * Delta(orb, n, spin1, spin2, lat, band)
              ! To get same critical temperature in all orbitals
              DO orb_prime = 1, discretization % ORBITALS
                IF (orb .NE. orb_prime) THEN
                  Gamma_SC(orb, n, spin1, spin2, lat, band) = Gamma_SC(orb, n, spin1, spin2, lat, band) - 0.5 * nearest_interorb_multiplier * Delta(orb_prime, n, spin1, spin2, lat, band)
                END IF
              END DO
              !To get same critical temperature in all subbands
              DO band_prime = 1, discretization % SUBBANDS
                IF (band .NE. band_prime) THEN
                  DO orb_prime = 1, discretization % ORBITALS
                    Gamma_SC(orb, n, spin1, spin2, lat, band) = Gamma_SC(orb, n, spin1, spin2, lat, band) - 0.5 * nearest_interorb_multiplier / 10.0d0 * Delta(orb_prime, n, spin1, spin2, lat, band_prime)
                  END DO
                END IF
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !Gamma for next nearest neighbours pairing
  DO band = 1, discretization % SUBBANDS
    DO spin1 = 1, SPINS   !Loop over spin coupling up-down or down-up
      DO spin2 = 1, SPINS
        DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS !Loop over next nearest neighbours
          DO lat = 1, discretization % SUBLATTICES !Loop over sublattices coupling, since for next-to-nearest neighbours it has intra-layer character
            DO orb = 1, discretization % ORBITALS

              Gamma_SC(orb, n, spin1, spin2, lat, band) = -0.5 * Delta(orb, n, spin1, spin2, lat, band)
              ! To get same critical temperature in all orbitals
              DO orb_prime = 1, discretization % ORBITALS
                IF (orb .NE. orb_prime) THEN
                  ! CHECK WHETHER THIS IS 0.25 OR 0.5
                  Gamma_SC(orb, n, spin1, spin2, lat, band) = Gamma_SC(orb, n, spin1, spin2, lat, band) - 0.5 * next_interorb_multiplier * Delta(orb_prime, n, spin1, spin2, lat, band)
                END IF
              END DO
              !To get same critical temperature in all subbands
              DO band_prime = 1, discretization % SUBBANDS
                IF (band .NE. band_prime) THEN
                  DO orb_prime = 1, discretization % ORBITALS
                    Gamma_SC(orb, n, spin1, spin2, lat, band) = Gamma_SC(orb, n, spin1, spin2, lat, band) - 0.5 * next_interorb_multiplier / 10.0d0 * Delta(orb_prime, n, spin1, spin2, lat, band_prime)
                  END DO
                END IF
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE GET_GAMMAS_FROM_DELTAS

SUBROUTINE CHECK_CONVERGENCE(sc_flag, Gamma_old, Gamma_new, Charge_dens, Charge_dens_new, gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error, sc_iter, self_consistency, discretization)
  IMPLICIT NONE

  TYPE(discretization_t), INTENT(IN) :: discretization
  TYPE(self_consistency_t), INTENT(INOUT) :: self_consistency
  LOGICAL, INTENT(INOUT) :: sc_flag
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(IN) :: Gamma_old(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(IN) :: Gamma_new(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#else
  COMPLEX(REAL64), INTENT(IN) :: Gamma_old(discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(IN) :: Gamma_new(discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#endif
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K, discretization % SUBBANDS)
  REAL(REAL64), INTENT(IN) :: Charge_dens_new(discretization % derived % DIM_POSITIVE_K, discretization % SUBBANDS)
  REAL(REAL64), INTENT(INOUT) :: gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error
  INTEGER(INT32), INTENT(IN) :: sc_iter

  gamma_max_error_prev = gamma_max_error
  charge_max_error_prev = charge_max_error
  gamma_max_error = 0.
  charge_max_error = 0.
  !Here we check whether convergence was reached
  sc_flag = .FALSE.
  gamma_max_error = MAXVAL(ABS(Gamma_new - Gamma_old))
  charge_max_error = MAXVAL(ABS(Charge_dens - Charge_dens_new))
  IF (gamma_max_error < self_consistency % gamma_eps_convergence .AND. charge_max_error < self_consistency % charge_eps_convergence) THEN
    sc_flag = .TRUE.
  END IF

  !Change Broyden mixing parameter if simulation diverges between iterations
  !To avoid oscillations near convergence
  IF (ABS(gamma_max_error) > ABS(gamma_max_error_prev) .AND. (sc_iter > 1)) THEN
    !PRINT*, "Adapted sc_alpha = ", sc_alpha
    WRITE (log_string, '(a, E15.5)') "Divergent iteration. Adapted sc_alpha: ", self_consistency % sc_alpha
    LOG_ABNORMAL(log_string)
    self_consistency % sc_alpha = self_consistency % sc_alpha * self_consistency % sc_alpha_adapt
  END IF
END SUBROUTINE CHECK_CONVERGENCE

SUBROUTINE FLATTEN_FOR_BROYDEN(Gamma_old, Gamma_new, Charge_dens, Charge_dens_new, Broyden_vector_new, Broyden_vector, broyden_length, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(IN) :: Gamma_old(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(IN) :: Gamma_new(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % SUBBANDS)
#else
  COMPLEX(REAL64), INTENT(IN) :: Gamma_old(discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(IN) :: Gamma_new(discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % derived % DIM_POSITIVE_K, &
                                          & discretization % SUBBANDS)
#endif
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K, discretization % SUBBANDS)
  REAL(REAL64), INTENT(IN) :: Charge_dens_new(discretization % derived % DIM_POSITIVE_K, discretization % SUBBANDS)
  INTEGER(INT32), INTENT(IN) :: broyden_length
  REAL(REAL64), INTENT(OUT) :: Broyden_vector_new(broyden_length)
  REAL(REAL64), INTENT(OUT) :: Broyden_vector(broyden_length)

  INTEGER(INT32) :: band, spin1, spin2, orb, n, lat, broyden_index
  INTEGER(INT32) :: broyden_index_new

  INTEGER(INT32) :: i_band, j_band, neigh
  broyden_index = 1
#ifndef BAND_BASIS
  DO band = 1, discretization % SUBBANDS
    DO i_band = 1, discretization % derived % DIM_POSITIVE_K
      DO j_band = 1, discretization % derived % DIM_POSITIVE_K
        DO neigh = 1, N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS
          Broyden_vector(broyden_index) = REAL(Gamma_old(neigh, i_band, j_band, band))
          Broyden_vector_new(broyden_index) = REAL(Gamma_new(neigh, i_band, j_band, band))
          broyden_index = broyden_index + 1

          Broyden_vector(broyden_index) = AIMAG(Gamma_old(neigh, i_band, j_band, band))
          Broyden_vector_new(broyden_index) = AIMAG(Gamma_new(neigh, i_band, j_band, band))
          broyden_index = broyden_index + 1
        END DO
      END DO
    END DO
  END DO
#else
  DO band = 1, discretization % SUBBANDS
    DO i_band = 1, discretization % derived % DIM_POSITIVE_K
      DO j_band = 1, discretization % derived % DIM_POSITIVE_K
        Broyden_vector(broyden_index) = REAL(Gamma_old(i_band, j_band, band))
        Broyden_vector_new(broyden_index) = REAL(Gamma_new(i_band, j_band, band))
        broyden_index = broyden_index + 1

        Broyden_vector(broyden_index) = AIMAG(Gamma_old(i_band, j_band, band))
        Broyden_vector_new(broyden_index) = AIMAG(Gamma_new(i_band, j_band, band))
        broyden_index = broyden_index + 1
      END DO
    END DO
  END DO

#endif

  !Must be +1!!!
  DO band = 1, discretization % SUBBANDS
    DO n = 1, discretization % derived % DIM_POSITIVE_K
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

SUBROUTINE RESHAPE_FROM_BROYDEN(Gamma_SC, Charge_dens, Broyden_vector, broyden_length, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(OUT) :: Gamma_SC(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#else
  COMPLEX(REAL64), INTENT(OUT) :: Gamma_SC(discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % derived % DIM_POSITIVE_K, &
                                         & discretization % SUBBANDS)
#endif
  REAL(REAL64), INTENT(OUT) :: Charge_dens(discretization % derived % DIM_POSITIVE_K, discretization % SUBBANDS)
  INTEGER(INT32), INTENT(IN) :: broyden_length
  REAL(REAL64), INTENT(IN) :: Broyden_vector(broyden_length)

  INTEGER(INT32) :: band, spin1, spin2, orb, n, lat, broyden_index
  INTEGER(INT32) :: i_band, j_band, neigh
  broyden_index = 1
#ifndef BAND_BASIS
  DO band = 1, discretization % SUBBANDS
    DO i_band = 1, discretization % derived % DIM_POSITIVE_K
      DO j_band = 1, discretization % derived % DIM_POSITIVE_K
        DO neigh = 1, N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS
          Gamma_SC(neigh, i_band, j_band, band) = CMPLX(Broyden_vector(broyden_index), Broyden_vector(broyden_index + 1), KIND=REAL64)
          broyden_index = broyden_index + 2
        END DO
      END DO
    END DO
  END DO
#else
  DO band = 1, discretization % SUBBANDS
    DO i_band = 1, discretization % derived % DIM_POSITIVE_K
      DO j_band = 1, discretization % derived % DIM_POSITIVE_K
        Gamma_SC(i_band, j_band, band) = CMPLX(Broyden_vector(broyden_index), Broyden_vector(broyden_index + 1), KIND=REAL64)
        broyden_index = broyden_index + 2
      END DO
    END DO
  END DO
#endif

  DO band = 1, discretization % SUBBANDS
    DO n = 1, discretization % derived % DIM_POSITIVE_K
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

END MODULE self_consistency
