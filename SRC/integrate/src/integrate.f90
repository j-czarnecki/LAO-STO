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

MODULE integrate
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
USE local_integrand
USE logger
USe types
IMPLICIT NONE

#include "macros_def.f90"

CONTAINS

!Adapted from "Numerical Recipes in Fortran Second Edition"
!William H. Press, Saul A. Teukolsky, W. T. Vetterling, B. P. Flannery

RECURSIVE SUBROUTINE ROMBERG_Y(Hamiltonian_const, Gamma_SC, Charge_dens, i_r, k1_steps, k2_chunk_min, k2_chunk_max, &
                    & Delta_local, Charge_dens_local, sc_input)

  IMPLICIT NONE
  TYPE(sc_input_params_t), INTENT(IN) :: sc_input
  COMPLEX(REAL64), INTENT(IN) :: Hamiltonian_const(sc_input % discretization % derived % DIM, &
                                             & sc_input % discretization % derived % DIM)
  INTEGER(INT32), INTENT(IN) :: i_r, k1_steps
  REAL(REAL64), INTENT(IN) :: k2_chunk_min, k2_chunk_max
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(sc_input % discretization % ORBITALS, &
                                    & N_ALL_NEIGHBOURS, &
                                    & SPINS, &
                                    & SPINS, &
                                    & sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(IN) :: Charge_dens(sc_input % discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64), INTENT(OUT) :: Delta_local(sc_input % discretization % ORBITALS, &
                                        & N_ALL_NEIGHBOURS, &
                                        & SPINS, &
                                        & SPINS, &
                                        & sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(OUT) :: Charge_dens_local(sc_input % discretization % derived % DIM_POSITIVE_K)
  !Parameters for Romberg integration

  COMPLEX(REAL64) :: stepsize(sc_input % romberg % max_grid_refinements_y + 1)
  COMPLEX(REAL64) :: Delta_iterations(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS, sc_input % romberg % max_grid_refinements_y + 1)
  REAL(REAL64) :: Charge_dens_iterations(sc_input % discretization % derived % DIM_POSITIVE_K, sc_input % romberg % max_grid_refinements_y + 1)
  COMPLEX(REAL64) :: Delta_sum(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64) :: Charge_dens_sum(sc_input % discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64) :: result_error, result
  INTEGER(INT32) :: n, i, j, spin1, spin2, orb, lat
  REAL(REAL64) :: dk2_trap, k2_trap
  LOGICAL :: convergence
  REAL(REAL64) :: r_max_local, k1_chunk_min, k1_chunk_max
  REAL(REAL64) :: max_error_delta, max_error_charge

  stepsize = CMPLX(0., 0., KIND=REAL64)
  Delta_iterations = CMPLX(0., 0., KIND=REAL64)
  Charge_dens_iterations = CMPLX(0., 0., KIND=REAL64)

  convergence = .FALSE.
  !stepsize(1) = k2_chunk_max - k2_chunk_min
  stepsize(1) = 1.
  DO j = 1, sc_input % romberg % max_grid_refinements_y
    !TRAPZD IMPLEMENTATION HERE
    !First approximation of the integral is taking only boundary values
    IF (j == 1) THEN
      !Calculation for lower bound of chunk
      r_max_local = r_max_phi(MOD(ABS(k2_chunk_max), PI / 3))
      k1_chunk_min = r_max_local / k1_steps * i_r
      k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
      CALL ROMBERG_X(Hamiltonian_const, Gamma_SC, Charge_dens, k1_chunk_min, k1_chunk_max, k2_chunk_max,&
          &  Delta_local, Charge_dens_local, sc_input)
      Delta_iterations(:, :, :, :, :, j) = Delta_local
      Charge_dens_iterations(:, j) = Charge_dens_local

      !Calculation for upper bound of chunk
      r_max_local = r_max_phi(MOD(ABS(k2_chunk_min), PI / 3))
      k1_chunk_min = r_max_local / k1_steps * i_r
      k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
      CALL ROMBERG_X(Hamiltonian_const, Gamma_SC, Charge_dens, k1_chunk_min, k1_chunk_max, k2_chunk_min,&
          &  Delta_local, Charge_dens_local, sc_input)
      Delta_iterations(:, :, :, :, :, j) = Delta_iterations(:, :, :, :, :, j) + Delta_local
      Charge_dens_iterations(:, j) = Charge_dens_iterations(:, j) + Charge_dens_local

      Delta_iterations(:, :, :, :, :, j) = 0.5 * (k2_chunk_max - k2_chunk_min) * Delta_iterations(:, :, :, :, :, j)
      Charge_dens_iterations(:, j) = 0.5 * (k2_chunk_max - k2_chunk_min) * Charge_dens_iterations(:, j)

      !Next approximations take point in between already calculated points
      ! i.e make the grid twice as dense as in previous iteration
    ELSE
      i = 2**(j - 2)
      dk2_trap = (k2_chunk_max - k2_chunk_min) / i
      k2_trap = k2_chunk_min + 0.5 * dk2_trap
      Delta_sum = CMPLX(0., 0., KIND=REAL64)
      Charge_dens_sum(:) = 0.
      ! Delta_iterations(:,:,:,:,j) = CMPLX(0. , 0., KIND=REAL64)
      ! Charge_dens_iterations(:,j) = CMPLX(0. , 0., KIND=REAL64)
      DO n = 1, i
        !Here we pass k1_trap as actual k1 point
        r_max_local = r_max_phi(MOD(ABS(k2_trap), PI / 3))
        k1_chunk_min = r_max_local / k1_steps * i_r
        k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
        CALL ROMBERG_X(Hamiltonian_const, Gamma_SC, Charge_dens, k1_chunk_min, k1_chunk_max, k2_trap,&
            &  Delta_local, Charge_dens_local, sc_input)
        Delta_sum = Delta_sum + Delta_local
        Charge_dens_sum = Charge_dens_sum + Charge_dens_local
        ! Delta_iterations(:,:,:,:,j) =  Delta_iterations(:,:,:,:,j) + Delta_local(:,:,:,:)
        ! Charge_dens_iterations(:,j) = Charge_dens_iterations(:,j) + Charge_dens_local(:)
        k2_trap = k2_trap + dk2_trap
      END DO
      Delta_iterations(:, :, :, :, :, j) = 0.5 * (Delta_iterations(:, :, :, :, :, j) + (k2_chunk_max - k2_chunk_min) * Delta_sum / i)
      Charge_dens_iterations(:, j) = 0.5 * (Charge_dens_iterations(:, j) + (k2_chunk_max - k2_chunk_min) * Charge_dens_sum / i)

      IF (j >= sc_input % romberg % interpolation_deg_y) THEN
        max_error_delta = 0.
        max_error_charge = 0.

        !For all components of Delta_iterations and Charge_dens_iterations
        !Check whether integral when dk1 ---> 0 can be approximated
        !With relative error no bigger than EPS
        convergence = .TRUE.
        !Checking Delta_iterations convergence
        DO spin1 = 1, SPINS
          DO spin2 = 1, SPINS
            DO orb = 1, sc_input % discretization % ORBITALS
              !Split into N_NEIGHBOURS and N_ALL_NEIGHBOURS because other parts of matrix are referenced
              DO n = 1, N_NEIGHBOURS
                DO lat = 1, sc_input % discretization % derived % LAYER_COUPLINGS
                  CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_y + 1):j), Delta_iterations(orb, n, spin1, spin2, lat, (j - sc_input % romberg % interpolation_deg_y + 1):j), &
                              & sc_input % romberg % interpolation_deg_y, CMPLX(0., 0., KIND=REAL64), result, result_error)
                  Delta_local(orb, n, spin1, spin2, lat) = result
                  IF (ABS(result) > 0.0d0) THEN
                    max_error_delta = MAX(max_error_delta, ABS(result_error) / ABS(result))
                  END IF
                  IF (ABS(result_error) > sc_input % romberg % romb_eps_y * ABS(result)) THEN
                    convergence = .FALSE.
                  END IF
                END DO
              END DO
              DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
                DO lat = 1, sc_input % discretization % SUBLATTICES
                  CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_y + 1):j), Delta_iterations(orb, n, spin1, spin2, lat, (j - sc_input % romberg % interpolation_deg_y + 1):j), &
                              & sc_input % romberg % interpolation_deg_y, CMPLX(0., 0., KIND=REAL64), result, result_error)
                  Delta_local(orb, n, spin1, spin2, lat) = result
                  IF (ABS(result) > 0.0d0) THEN
                    max_error_delta = MAX(max_error_delta, ABS(result_error) / ABS(result))
                  END IF
                  IF (ABS(result_error) > sc_input % romberg % romb_eps_y * ABS(result)) THEN
                    convergence = .FALSE.
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO

        !Checking Charge_dens convergence
        DO n = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_y + 1):j), CMPLX(Charge_dens_iterations(n, (j - sc_input % romberg % interpolation_deg_y + 1):j), 0.0_REAL64, KIND=REAL64), &
                      & sc_input % romberg % interpolation_deg_y, CMPLX(0., 0., KIND=REAL64), result, result_error)
          Charge_dens_local(n) = REAL(result)
          max_error_charge = MAX(max_error_charge, ABS(result_error) / ABS(result))
          IF (ABS(result_error) > sc_input % romberg % romb_eps_y * ABS(result)) THEN
            convergence = .FALSE.
          END IF
        END DO

      END IF
    END IF

    IF (convergence) THEN
      WRITE (log_string, '(a, I15)') "Romberg Y converged after iteration: ", j
      LOG_DEBUG(log_string)
      RETURN
    ELSE
      Delta_iterations(:, :, :, :, :, j + 1) = Delta_iterations(:, :, :, :, :, j)
      Charge_dens_iterations(:, j + 1) = Charge_dens_iterations(:, j)
      stepsize(j + 1) = 0.25 * stepsize(j)
    END IF
  END DO

  WRITE (log_string, '(2(a, F10.6), a, I15, 2(a, E15.5))') "Romberg Y did not converge for chunk &
  & k1_chunk_min: ", k1_chunk_min, " k2_chunk_min: ", k2_chunk_min, " after iteration: ", j - 1, &
  & " with max_error_delta: ", max_error_delta, " max_error_charge: ", max_error_charge
  LOG_ABNORMAL(log_string)

END SUBROUTINE ROMBERG_Y

RECURSIVE SUBROUTINE ROMBERG_X(Hamiltonian_const, Gamma_SC, Charge_dens, k1_chunk_min, k1_chunk_max, k2_actual, Delta_local, Charge_dens_local, sc_input)
  IMPLICIT NONE
  TYPE(sc_input_params_t), INTENT(IN) :: sc_input
  COMPLEX(REAL64), INTENT(IN) :: Hamiltonian_const(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: k1_chunk_min, k1_chunk_max, k2_actual
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(IN) :: Charge_dens(sc_input % discretization % derived % DIM_POSITIVE_K)

  COMPLEX(REAL64), INTENT(OUT) :: Delta_local(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(OUT) :: Charge_dens_local(sc_input % discretization % derived % DIM_POSITIVE_K)

  COMPLEX(REAL64) :: stepsize(sc_input % romberg % max_grid_refinements_x + 1)
  COMPLEX(REAL64) :: Delta_iterations(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS, sc_input % romberg % max_grid_refinements_x + 1)
  REAL(REAL64) :: Charge_dens_iterations(sc_input % discretization % derived % DIM_POSITIVE_K, sc_input % romberg % max_grid_refinements_x + 1)
  COMPLEX(REAL64) :: Delta_sum(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64) :: Charge_dens_sum(sc_input % discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64) :: result_error, result
  INTEGER(INT32) :: n, i, j, spin1, spin2, orb, lat
  REAL(REAL64) :: dk1_trap, k1_trap
  LOGICAL :: convergence

  stepsize = CMPLX(0., 0., KIND=REAL64)
  Delta_iterations = CMPLX(0., 0., KIND=REAL64)
  Charge_dens_iterations = CMPLX(0., 0., KIND=REAL64)

  convergence = .FALSE.
  !stepsize(1) = k1_chunk_max - k1_chunk_min
  stepsize(1) = 1.

  DO j = 1, sc_input % romberg % max_grid_refinements_x
    !TRAPZD IMPLEMENTATION HERE
    !First approximation of the integral is taking only boundary values
    IF (j == 1) THEN
      !Calculation for lower bound of chunk
      CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, &
          & Charge_dens, k1_chunk_min, k2_actual, Delta_local, Charge_dens_local, sc_input % discretization, sc_input % physical)
      Delta_iterations(:, :, :, :, :, j) = Delta_local
      Charge_dens_iterations(:, j) = Charge_dens_local(:)

      !Calculation for upper bound of chunk
      CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, &
      & Charge_dens, k1_chunk_max, k2_actual, Delta_local, Charge_dens_local, sc_input % discretization, sc_input % physical)
      Delta_iterations(:, :, :, :, :, j) = Delta_iterations(:, :, :, :, :, j) + Delta_local
      Charge_dens_iterations(:, j) = Charge_dens_iterations(:, j) + Charge_dens_local

      Delta_iterations(:, :, :, :, :, j) = 0.5 * (k1_chunk_max - k1_chunk_min) * Delta_iterations(:, :, :, :, :, j)
      Charge_dens_iterations(:, j) = 0.5 * (k1_chunk_max - k1_chunk_min) * Charge_dens_iterations(:, j)

      !Next approximations take point in between already calculated points
      ! i.e make the grid twice as dense as in previous iteration
    ELSE
      i = 2**(j - 2)
      dk1_trap = (k1_chunk_max - k1_chunk_min) / i
      k1_trap = k1_chunk_min + 0.5 * dk1_trap
      Delta_sum = CMPLX(0., 0., KIND=REAL64)
      Charge_dens_sum = 0.
      DO n = 1, i
        !Here we pass k1_trap as actual k1 point
        CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, &
        & Charge_dens, k1_trap, k2_actual, Delta_local, Charge_dens_local, sc_input % discretization, sc_input % physical)
        Delta_sum = Delta_sum + Delta_local
        Charge_dens_sum = Charge_dens_sum + Charge_dens_local
        k1_trap = k1_trap + dk1_trap
      END DO

      Delta_iterations(:, :, :, :, :, j) = 0.5 * (Delta_iterations(:, :, :, :, :, j) + (k1_chunk_max - k1_chunk_min) * Delta_sum / i)
      Charge_dens_iterations(:, j) = 0.5 * (Charge_dens_iterations(:, j) + (k1_chunk_max - k1_chunk_min) * Charge_dens_sum / i)

      IF (j >= sc_input % romberg % interpolation_deg_x) THEN
        !For all components of Delta_iterations and Charge_dens_iterations
        !Check whether integral when dk1 ---> 0 can be approximated
        !With relative error no bigger than EPS
        convergence = .TRUE.
        !Checking Delta_iterations convergence
        DO spin1 = 1, SPINS
          DO spin2 = 1, SPINS
            DO orb = 1, sc_input % discretization % ORBITALS
              !Split into N_NEIGHBOURS and N_ALL_NEIGHBOURS because other parts of matrix are referenced
              DO n = 1, N_NEIGHBOURS
                DO lat = 1, sc_input % discretization % derived % LAYER_COUPLINGS
                  CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_x + 1):j), Delta_iterations(orb, n, spin1, spin2, lat, (j - sc_input % romberg % interpolation_deg_x + 1):j), &
                              & sc_input % romberg % interpolation_deg_x, CMPLX(0., 0., KIND=REAL64), result, result_error)
                  Delta_local(orb, n, spin1, spin2, lat) = result
                  IF (ABS(result_error) > sc_input % romberg % romb_eps_x * ABS(result)) THEN
                    convergence = .FALSE.
                  END IF
                END DO
              END DO
              DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
                DO lat = 1, sc_input % discretization % SUBLATTICES
                  CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_x + 1):j), Delta_iterations(orb, n, spin1, spin2, lat, (j - sc_input % romberg % interpolation_deg_x + 1):j), &
                              & sc_input % romberg % interpolation_deg_x, CMPLX(0., 0., KIND=REAL64), result, result_error)
                  Delta_local(orb, n, spin1, spin2, lat) = result
                  IF (ABS(result_error) > sc_input % romberg % romb_eps_x * ABS(result)) THEN
                    convergence = .FALSE.
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO
        !Checking Charge_dens convergence
        DO n = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          CALL POLINT(stepsize((j - sc_input % romberg % interpolation_deg_x + 1):j), CMPLX(Charge_dens_iterations(n, (j - sc_input % romberg % interpolation_deg_x + 1):j), 0.0_REAL64, KIND=REAL64), &
                      & sc_input % romberg % interpolation_deg_x, CMPLX(0., 0., KIND=REAL64), result, result_error)
          Charge_dens_local(n) = REAL(result)
          IF (ABS(result_error) > sc_input % romberg % romb_eps_x * ABS(result)) THEN
            convergence = .FALSE.
          END IF
        END DO
      END IF
    END IF

    IF (convergence) THEN
      !No log necessary, since Romberg Y will give info about this
      WRITE (log_string, '(a, I15)') "Romberg X converged after iteration: ", j
      LOG_DEBUG(log_string)
      RETURN
    ELSE
      Delta_iterations(:, :, :, :, :, j + 1) = Delta_iterations(:, :, :, :, :, j)
      Charge_dens_iterations(:, j + 1) = Charge_dens_iterations(:, j)
      stepsize(j + 1) = 0.25 * stepsize(j)
    END IF

  END DO

  ! WRITE(log_string,'(a, F10.6, a, F10.6, a, I15)') "Romberg X did not converge for &
  ! & k1_chunk_min: ", k1_chunk_min, " k2_actual: ", k2_actual, " after iteration: ", j - 1
  ! LOG_ABNORMAL(log_string)

END SUBROUTINE ROMBERG_X

!See Section 3.1
RECURSIVE SUBROUTINE POLINT(X, Y, deg, x_target, y_approx, dy)
  INTEGER(INT32), INTENT(IN) :: deg
  COMPLEX(REAL64), INTENT(IN) :: X(deg), Y(deg), x_target
  COMPLEX(REAL64), INTENT(OUT) :: y_approx, dy

  INTEGER(INT32) :: i, m, nearest
  COMPLEX(REAL64) :: dx_left, dx_right, w, den
  REAL(REAL64) :: diff, diff_temp
  COMPLEX(REAL64) :: C(deg), D(deg)

  nearest = 1
  diff = ABS(x_target - X(1))

  !Finding tabulated X closest to x_target
  DO i = 1, nearest
    diff_temp = ABS(x_target - X(i))
    IF (diff_temp < diff) THEN
      diff = diff_temp
      nearest = i
    END IF
  END DO

  !Initializing C and D values
  C(:) = Y(:)
  D(:) = Y(:)
  y_approx = Y(nearest)
  nearest = nearest - 1

  DO m = 1, deg - 1
    DO i = 1, deg - m
      dx_left = X(i) - x_target
      dx_right = X(i + m) - x_target
      w = C(i + 1) - D(i)
      den = dx_left - dx_right
      IF (den == 0) STOP 'Repeating values in X table'
      den = w / den
      D(i) = dx_right * den
      C(i) = dx_left * den
    END DO

    IF (2 * nearest < deg - m) THEN
      dy = C(nearest + 1)
    ELSE
      dy = D(nearest)
      nearest = nearest - 1
    END IF
    y_approx = y_approx + dy
  END DO

END SUBROUTINE

! !See section 4.3
! !Romberg integration
! SUBROUTINE QROMB(func, x_min, x_max, result)
!     COMPLEX(REAL64), INTENT(IN) :: x_min, x_max
!     COMPLEX(REAL64), EXTERNAL :: func
!     COMPLEX(REAL64), INTENT(OUT) :: result

!     REAL(REAL64), PARAMETER :: EPS = 1e-6
!     INTEGER(INT32), PARAMETER :: JMAX = 20
!     INTEGER(INT32), PARAMETER :: K = 5

!     COMPLEX(REAL64) :: stepsize(JMAX + 1), integral(JMAX + 1)
!     COMPLEX(REAL64) :: result_error
!     INTEGER(INT32) :: j

!     stepsize(1) = 1.
!     !Maximum number of splittings to half
!     DO j = 1, JMAX
!         CALL TRAPZD(func, x_min, x_max, integral(j), j)
!         IF (j >= K) THEN
!             !Really smart - based on previous iterations we want to extrapolate value of integral when stepsize would be 0.
!             !Passing last K steps
!             CALL POLINT(stepsize(j - K + 1), integral(j - K + 1), K, CMPLX(0. , 0., KIND=REAL64), result, result_error)
!             IF (ABS(result_error) < EPS*ABS(result)) RETURN
!         END IF
!         integral(j + 1) = integral(j)
!         stepsize(j+1) = 0.25*stepsize(j)    !This is crucial
!     END DO

! END SUBROUTINE QROMB

! !See Section 4.2
! !Extended trapezoidal rule
! SUBROUTINE TRAPZD(func, x_min, x_max, result, n)
!     COMPLEX(REAL64), INTENT(IN) :: x_max, x_min
!     INTEGER(INT32), INTENT(IN) :: n
!     COMPLEX(REAL64), INTENT(OUT) :: result
!     COMPLEX(REAL64), EXTERNAL :: func

!     INTEGER(INT32) :: i,j
!     COMPLEX(REAL64) :: dx, x, sum

!     !First approximation of integral
!     IF (n == 1) THEN
!         result = 0.5*(x_max - x_min)*(func(x_min) + func(x_max))
!     ELSE
!         i = 2**(n-2)
!         dx = (x_max  - x_min)/i
!         x = x_min + 0.5*dx
!         sum = 0.
!         DO j = 1, i
!             sum = sum + func(x)
!             x = x + dx
!         END DO
!         result = 0.5*(result + (x_max - x_min)*sum/i)
!     END IF

! END SUBROUTINE TRAPZD

END MODULE integrate
