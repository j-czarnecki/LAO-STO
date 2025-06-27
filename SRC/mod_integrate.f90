MODULE mod_integrate
USE mod_parameters
USE mod_local_integrand
USE mod_logger
IMPLICIT NONE

#include "macros_def.f90"

CONTAINS

!Adapted from "Numerical Recipes in Fortran Second Edition"
!William H. Press, Saul A. Teukolsky, W. T. Vetterling, B. P. Flannery

RECURSIVE SUBROUTINE ROMBERG_Y(Hamiltonian_const, Gamma_SC, Charge_dens, i_r, k1_steps, k2_chunk_min, k2_chunk_max, &
                    & Delta_local, Charge_dens_local, romb_eps_x, interpolation_deg_x, max_grid_refinements_x, &
                    & romb_eps_y, interpolation_deg_y, max_grid_refinements_y)

  COMPLEX*16, INTENT(IN) :: Hamiltonian_const(DIM, DIM)
  INTEGER*4, INTENT(IN) :: i_r, k1_steps
  REAL*8, INTENT(IN) :: k2_chunk_min, k2_chunk_max
  COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K)
  REAL*8, INTENT(IN) :: romb_eps_x, romb_eps_y
  INTEGER*4, INTENT(IN) :: interpolation_deg_x, interpolation_deg_y, max_grid_refinements_x, max_grid_refinements_y
  COMPLEX*16, INTENT(OUT) :: Delta_local(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(OUT) :: Charge_dens_local(DIM_POSITIVE_K)
  !Parameters for Romberg integration

  COMPLEX*16 :: stepsize(max_grid_refinements_y + 1)
  COMPLEX*16 :: Delta_iterations(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, max_grid_refinements_y + 1)
  REAL*8 :: Charge_dens_iterations(DIM_POSITIVE_K, max_grid_refinements_y + 1)
  COMPLEX*16 :: Delta_sum(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8 :: Charge_dens_sum(DIM_POSITIVE_K)
  COMPLEX*16 :: result_error, result
  INTEGER*4 :: n, i, j, spin, orb, lat
  REAL*8 :: dk2_trap, k2_trap
  LOGICAL :: convergence
  REAL*8 :: r_max_local, k1_chunk_min, k1_chunk_max

  stepsize = DCMPLX(0., 0.)
  Delta_iterations = DCMPLX(0., 0.)
  Charge_dens_iterations = DCMPLX(0., 0.)

  convergence = .FALSE.
  !stepsize(1) = k2_chunk_max - k2_chunk_min
  stepsize(1) = 1.
  DO j = 1, max_grid_refinements_y
    !TRAPZD IMPLEMENTATION HERE
    !First approximation of the integral is taking only boundary values
    IF (j == 1) THEN
      !Calculation for lower bound of chunk
      r_max_local = r_max_phi(MOD(ABS(k2_chunk_max), PI / 3))
      k1_chunk_min = r_max_local / k1_steps * i_r
      k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
      CALL ROMBERG_X(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), Charge_dens(:), k1_chunk_min, k1_chunk_max, k2_chunk_max,&
          &  Delta_local(:, :, :, :), Charge_dens_local(:), romb_eps_x, interpolation_deg_x, max_grid_refinements_x)
      Delta_iterations(:, :, :, :, j) = Delta_local(:, :, :, :)
      Charge_dens_iterations(:, j) = Charge_dens_local(:)

      !Calculation for upper bound of chunk
      r_max_local = r_max_phi(MOD(ABS(k2_chunk_min), PI / 3))
      k1_chunk_min = r_max_local / k1_steps * i_r
      k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
      CALL ROMBERG_X(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), Charge_dens(:), k1_chunk_min, k1_chunk_max, k2_chunk_min,&
          &  Delta_local(:, :, :, :), Charge_dens_local(:), romb_eps_x, interpolation_deg_x, max_grid_refinements_x)
      Delta_iterations(:, :, :, :, j) = Delta_iterations(:, :, :, :, j) + Delta_local(:, :, :, :)
      Charge_dens_iterations(:, j) = Charge_dens_iterations(:, j) + Charge_dens_local(:)

      Delta_iterations(:, :, :, :, j) = 0.5 * (k2_chunk_max - k2_chunk_min) * Delta_iterations(:, :, :, :, j)
      Charge_dens_iterations(:, j) = 0.5 * (k2_chunk_max - k2_chunk_min) * Charge_dens_iterations(:, j)

      !Next approximations take point in between already calculated points
      ! i.e make the grid twice as dense as in previous iteration
    ELSE
      i = 2**(j - 2)
      dk2_trap = (k2_chunk_max - k2_chunk_min) / i
      k2_trap = k2_chunk_min + 0.5 * dk2_trap
      Delta_sum(:, :, :, :) = DCMPLX(0., 0.)
      Charge_dens_sum(:) = 0.
      ! Delta_iterations(:,:,:,:,j) = DCMPLX(0. , 0.)
      ! Charge_dens_iterations(:,j) = DCMPLX(0. , 0.)
      DO n = 1, i
        !Here we pass k1_trap as actual k1 point
        r_max_local = r_max_phi(MOD(ABS(k2_trap), PI / 3))
        k1_chunk_min = r_max_local / k1_steps * i_r
        k1_chunk_max = r_max_local / k1_steps * (i_r + 1)
        CALL ROMBERG_X(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), Charge_dens(:), k1_chunk_min, k1_chunk_max, k2_trap,&
            &  Delta_local(:, :, :, :), Charge_dens_local(:), romb_eps_x, interpolation_deg_x, max_grid_refinements_x)
        Delta_sum(:, :, :, :) = Delta_sum(:, :, :, :) + Delta_local(:, :, :, :)
        Charge_dens_sum(:) = Charge_dens_sum(:) + Charge_dens_local(:)
        ! Delta_iterations(:,:,:,:,j) =  Delta_iterations(:,:,:,:,j) + Delta_local(:,:,:,:)
        ! Charge_dens_iterations(:,j) = Charge_dens_iterations(:,j) + Charge_dens_local(:)
        k2_trap = k2_trap + dk2_trap
      END DO
      Delta_iterations(:, :, :, :, j) = 0.5 * (Delta_iterations(:, :, :, :, j) + (k2_chunk_max - k2_chunk_min) * Delta_sum(:, :, :, :) / i)
      Charge_dens_iterations(:, j) = 0.5 * (Charge_dens_iterations(:, j) + (k2_chunk_max - k2_chunk_min) * Charge_dens_sum(:) / i)

            !!!!!!!! THIS PART COULD LEAD TO A PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Delta_iterations(:,:,:,:,j) = 0.5*(k2_chunk_max - k2_chunk_min)/i * Delta_iterations(:,:,:,:,j)
      ! Charge_dens_iterations(:,j) = 0.5*(k2_chunk_max - k2_chunk_min)/i * Charge_dens_iterations(:,j)
      !Summing all previously calculated values
      !They were calculated in different points
      ! DO n = 1, j
      !     Delta_iterations(:,:,:,:,j) =  Delta_iterations(:,:,:,:,j) + 0.5*Delta_iterations(:,:,:,:,n)
      !     Charge_dens_iterations(:,j) = Charge_dens_iterations(:,j) + 0.5*Charge_dens_iterations(:,n)
      ! END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (j >= interpolation_deg_y) THEN
        !For all components of Delta_iterations and Charge_dens_iterations
        !Check whether integral when dk1 ---> 0 can be approximated
        !With relative error no bigger than EPS
        convergence = .TRUE.
        !Checking Delta_iterations convergence
        DO spin = 1, 2
          DO orb = 1, ORBITALS
            !Split into N_NEIGHBOURS and N_ALL_NEIGHBOURS because other parts of matrix are referenced
            DO n = 1, N_NEIGHBOURS
              DO lat = 1, LAYER_COUPLINGS
                CALL POLINT(stepsize((j - interpolation_deg_y + 1):j), Delta_iterations(orb, n, spin, lat, (j - interpolation_deg_y + 1):j), &
                            & interpolation_deg_y, DCMPLX(0., 0.), result, result_error)
                Delta_local(orb, n, spin, lat) = result
                IF (ABS(result_error) > romb_eps_y * ABS(result)) THEN
                  convergence = .FALSE.
                END IF
              END DO
            END DO
            DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
              DO lat = 1, SUBLATTICES
                CALL POLINT(stepsize((j - interpolation_deg_y + 1):j), Delta_iterations(orb, n, spin, lat, (j - interpolation_deg_y + 1):j), &
                            & interpolation_deg_y, DCMPLX(0., 0.), result, result_error)
                Delta_local(orb, n, spin, lat) = result
                IF (ABS(result_error) > romb_eps_y * ABS(result)) THEN
                  convergence = .FALSE.
                END IF
              END DO
            END DO
          END DO
        END DO

        !Chaecking Charge_dens convergence
        DO n = 1, DIM_POSITIVE_K
          CALL POLINT(stepsize((j - interpolation_deg_y + 1):j), DCMPLX(Charge_dens_iterations(n, (j - interpolation_deg_y + 1):j), 0.), &
                      & interpolation_deg_y, DCMPLX(0., 0.), result, result_error)
          Charge_dens_local(n) = REAL(result)
          IF (ABS(result_error) > romb_eps_y * ABS(result)) THEN
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
      Delta_iterations(:, :, :, :, j + 1) = Delta_iterations(:, :, :, :, j)
      Charge_dens_iterations(:, j + 1) = Charge_dens_iterations(:, j)
      stepsize(j + 1) = 0.25 * stepsize(j)
    END IF
  END DO

  WRITE (log_string, '(a, F10.6, a, F10.6, a, I15)') "Romberg Y did not converge for chunk &
  & k1_chunk_min: ", k1_chunk_min, " k2_chunk_min: ", k2_chunk_min, " after iteration: ", j - 1
  LOG_ABNORMAL(log_string)

END SUBROUTINE ROMBERG_Y

RECURSIVE SUBROUTINE ROMBERG_X(Hamiltonian_const, Gamma_SC, Charge_dens, k1_chunk_min, k1_chunk_max, k2_actual, Delta_local, Charge_dens_local, &
                    & romb_eps_x, interpolation_deg_x, max_grid_refinements_x)
  COMPLEX*16, INTENT(IN) :: Hamiltonian_const(DIM, DIM)
  REAL*8, INTENT(IN) :: k1_chunk_min, k1_chunk_max, k2_actual
  COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K)
  REAL*8, INTENT(IN) :: romb_eps_x
  INTEGER*4, INTENT(IN) :: interpolation_deg_x, max_grid_refinements_x

  COMPLEX*16, INTENT(OUT) :: Delta_local(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(OUT) :: Charge_dens_local(DIM_POSITIVE_K)

  COMPLEX*16 :: stepsize(max_grid_refinements_x + 1)
  COMPLEX*16 :: Delta_iterations(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, max_grid_refinements_x + 1)
  REAL*8 :: Charge_dens_iterations(DIM_POSITIVE_K, max_grid_refinements_x + 1)
  COMPLEX*16 :: Delta_sum(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8 :: Charge_dens_sum(DIM_POSITIVE_K)
  COMPLEX*16 :: result_error, result
  INTEGER*4 :: n, i, j, spin, orb, lat
  REAL*8 :: dk1_trap, k1_trap
  LOGICAL :: convergence

  stepsize = DCMPLX(0., 0.)
  Delta_iterations = DCMPLX(0., 0.)
  Charge_dens_iterations = DCMPLX(0., 0.)

  convergence = .FALSE.
  !stepsize(1) = k1_chunk_max - k1_chunk_min
  stepsize(1) = 1.

  DO j = 1, max_grid_refinements_x
    !TRAPZD IMPLEMENTATION HERE
    !First approximation of the integral is taking only boundary values
    IF (j == 1) THEN
      !Calculation for lower bound of chunk
      CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), &
          & Charge_dens(:), k1_chunk_min, k2_actual, Delta_local, Charge_dens_local)
      Delta_iterations(:, :, :, :, j) = Delta_local(:, :, :, :)
      Charge_dens_iterations(:, j) = Charge_dens_local(:)

      !Calculation for upper bound of chunk
      CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), &
      & Charge_dens(:), k1_chunk_max, k2_actual, Delta_local, Charge_dens_local)
      Delta_iterations(:, :, :, :, j) = Delta_iterations(:, :, :, :, j) + Delta_local(:, :, :, :)
      Charge_dens_iterations(:, j) = Charge_dens_iterations(:, j) + Charge_dens_local(:)

      Delta_iterations(:, :, :, :, j) = 0.5 * (k1_chunk_max - k1_chunk_min) * Delta_iterations(:, :, :, :, j)
      Charge_dens_iterations(:, j) = 0.5 * (k1_chunk_max - k1_chunk_min) * Charge_dens_iterations(:, j)

      !Next approximations take point in between already calculated points
      ! i.e make the grid twice as dense as in previous iteration
    ELSE
      i = 2**(j - 2)
      dk1_trap = (k1_chunk_max - k1_chunk_min) / i
      k1_trap = k1_chunk_min + 0.5 * dk1_trap
      Delta_sum(:, :, :, :) = DCMPLX(0., 0.)
      Charge_dens_sum(:) = 0.
      DO n = 1, i
        !Here we pass k1_trap as actual k1 point
        CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const(:, :), Gamma_SC(:, :, :, :), &
        & Charge_dens(:), k1_trap, k2_actual, Delta_local, Charge_dens_local)
        Delta_sum(:, :, :, :) = Delta_sum(:, :, :, :) + Delta_local(:, :, :, :)
        Charge_dens_sum(:) = Charge_dens_sum(:) + Charge_dens_local(:)
        ! Delta_iterations(:,:,:,:,j) =  Delta_iterations(:,:,:,:,j) + Delta_local(:,:,:,:)
        ! Charge_dens_iterations(:,j) = Charge_dens_iterations(:,j) + Charge_dens_local(:)
        k1_trap = k1_trap + dk1_trap
      END DO

      Delta_iterations(:, :, :, :, j) = 0.5 * (Delta_iterations(:, :, :, :, j) + (k1_chunk_max - k1_chunk_min) * Delta_sum(:, :, :, :) / i)
      Charge_dens_iterations(:, j) = 0.5 * (Charge_dens_iterations(:, j) + (k1_chunk_max - k1_chunk_min) * Charge_dens_sum(:) / i)

            !!!!!!!! THIS PART COULD LEAD TO A PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Delta_iterations(:,:,:,:,j) = 0.5*(k1_chunk_max - k1_chunk_min)/i * Delta_iterations(:,:,:,:,j)
      ! Charge_dens_iterations(:,j) = 0.5*(k1_chunk_max - k1_chunk_min)/i * Charge_dens_iterations(:,j)
      !Summing all previously calculated values
      !They were calculated in different points
      ! DO n = 1, j
      !     Delta_iterations(:,:,:,:,j) =  Delta_iterations(:,:,:,:,j) + 0.5*Delta_iterations(:,:,:,:,n)
      !     Charge_dens_iterations(:,j) = Charge_dens_iterations(:,j) + 0.5*Charge_dens_iterations(:,n)
      ! END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (j >= interpolation_deg_x) THEN
        !For all components of Delta_iterations and Charge_dens_iterations
        !Check whether integral when dk1 ---> 0 can be approximated
        !With relative error no bigger than EPS
        convergence = .TRUE.
        !Checking Delta_iterations convergence
        DO spin = 1, 2
          DO orb = 1, ORBITALS
            !Split into N_NEIGHBOURS and N_ALL_NEIGHBOURS because other parts of matrix are referenced
            DO n = 1, N_NEIGHBOURS
              DO lat = 1, LAYER_COUPLINGS
                CALL POLINT(stepsize((j - interpolation_deg_x + 1):j), Delta_iterations(orb, n, spin, lat, (j - interpolation_deg_x + 1):j), &
                            & interpolation_deg_x, DCMPLX(0., 0.), result, result_error)
                Delta_local(orb, n, spin, lat) = result
                IF (ABS(result_error) > romb_eps_x * ABS(result)) THEN
                  convergence = .FALSE.
                END IF
              END DO
            END DO
            DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
              DO lat = 1, SUBLATTICES
                CALL POLINT(stepsize((j - interpolation_deg_x + 1):j), Delta_iterations(orb, n, spin, lat, (j - interpolation_deg_x + 1):j), &
                            & interpolation_deg_x, DCMPLX(0., 0.), result, result_error)
                Delta_local(orb, n, spin, lat) = result
                IF (ABS(result_error) > romb_eps_x * ABS(result)) THEN
                  convergence = .FALSE.
                END IF
              END DO
            END DO
          END DO
        END DO

        !Chaecking Charge_dens convergence
        DO n = 1, DIM_POSITIVE_K
          CALL POLINT(stepsize((j - interpolation_deg_x + 1):j), DCMPLX(Charge_dens_iterations(n, (j - interpolation_deg_x + 1):j), 0.), &
                      & interpolation_deg_x, DCMPLX(0., 0.), result, result_error)
          Charge_dens_local(n) = REAL(result)
          IF (ABS(result_error) > romb_eps_x * ABS(result)) THEN
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
      Delta_iterations(:, :, :, :, j + 1) = Delta_iterations(:, :, :, :, j)
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
  INTEGER*4, INTENT(IN) :: deg
  COMPLEX*16, INTENT(IN) :: X(deg), Y(deg), x_target
  COMPLEX*16, INTENT(OUT) :: y_approx, dy

  INTEGER*4 :: i, m, nearest
  COMPLEX*16 :: dx_left, dx_right, w, den
  REAL*8 :: diff, diff_temp
  COMPLEX*16 :: C(deg), D(deg)

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
!     COMPLEX*16, INTENT(IN) :: x_min, x_max
!     COMPLEX*16, EXTERNAL :: func
!     COMPLEX*16, INTENT(OUT) :: result

!     REAL*8, PARAMETER :: EPS = 1e-6
!     INTEGER*4, PARAMETER :: JMAX = 20
!     INTEGER*4, PARAMETER :: K = 5

!     COMPLEX*16 :: stepsize(JMAX + 1), integral(JMAX + 1)
!     COMPLEX*16 :: result_error
!     INTEGER*4 :: j

!     stepsize(1) = 1.
!     !Maximum number of splittings to half
!     DO j = 1, JMAX
!         CALL TRAPZD(func, x_min, x_max, integral(j), j)
!         IF (j >= K) THEN
!             !Really smart - based on previous iterations we want to extrapolate value of integral when stepsize would be 0.
!             !Passing last K steps
!             CALL POLINT(stepsize(j - K + 1), integral(j - K + 1), K, DCMPLX(0. , 0.), result, result_error)
!             IF (ABS(result_error) < EPS*ABS(result)) RETURN
!         END IF
!         integral(j + 1) = integral(j)
!         stepsize(j+1) = 0.25*stepsize(j)    !This is crucial
!     END DO

! END SUBROUTINE QROMB

! !See Section 4.2
! !Extended trapezoidal rule
! SUBROUTINE TRAPZD(func, x_min, x_max, result, n)
!     COMPLEX*16, INTENT(IN) :: x_max, x_min
!     INTEGER*4, INTENT(IN) :: n
!     COMPLEX*16, INTENT(OUT) :: result
!     COMPLEX*16, EXTERNAL :: func

!     INTEGER*4 :: i,j
!     COMPLEX*16 :: dx, x, sum

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

END MODULE mod_integrate
