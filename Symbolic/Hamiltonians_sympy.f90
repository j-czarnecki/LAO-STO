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

MODULE hamiltonians_sympy
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USe omp_lib
IMPLICIT NONE

PRIVATE
PUBLIC :: COMPUTE_DOS, COMPUTE_FREE_ENERGY

INTEGER(INT32), PARAMETER :: DIM_POSITIVE_K = 12
REAL(REAL64), PARAMETER :: PI = 4 * ATAN(1.0d0)
! Integration constants
REAL(REAL64), PARAMETER :: K1_MAX = 1. !Full Brillouin zone to integrate over
REAL(REAL64), PARAMETER :: K2_MAX = 1.
REAL(REAL64), PARAMETER :: KX_MAX = 4.0d0 * PI / (3.0d0 * SQRT(3.0d0)) ! This defines maximum kx of hexagon
! corresponding to the first Brillouin Zone
REAL(REAL64), PARAMETER :: KY_MAX = 2.0d0 * PI / 3.0d0 ! This defines maximum ky of hexagon
! corresponding to the first Brillouin Zone
REAL(REAL64), PARAMETER :: R_K_MAX = 4.0d0 * PI / (3.0d0 * SQRT(3.0d0)) ! This defines the radius of circle
! that the first Brillouin Zone hexagon is inscribed in.
REAL(REAL64), PARAMETER :: JACOBIAN = 8 * PI**2 / (3.*SQRT(3.0d0))
INTEGER(INT32), PARAMETER :: N_BZ_SECTIONS = 6
REAL(REAL64), PARAMETER :: k_B = 8.617333262 * 1e-5 * 1e3 ! [meV]

CONTAINS

RECURSIVE SUBROUTINE DIAGONALIZE_HERMITIAN(Hamiltonian, Eigenvalues, N)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: N
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(N, N)
  REAL(REAL64), INTENT(OUT) :: Eigenvalues(N)
  COMPLEX(REAL64), ALLOCATABLE :: WORK(:)
  REAL(REAL64), ALLOCATABLE :: RWORK(:)
  INTEGER(INT32) :: LWORK
  INTEGER(INT32) :: INFO
  LWORK = 10 * N
  ALLOCATE (WORK(LWORK))
  ALLOCATE (RWORK(3 * N - 2))
  CALL ZHEEV('V', 'U', N, Hamiltonian, N, Eigenvalues, WORK, LWORK, RWORK, INFO)
  !WRITE(*, fmt = "(15F10.4)") Eigenvalues
  IF (INFO .ne. 0) THEN
    PRINT *, 'ZHEEV INFO ', INFO
    STOP
  END IF

  DEALLOCATE (WORK)
  DEALLOCATE (RWORK)

END SUBROUTINE DIAGONALIZE_HERMITIAN

!dir$ attributes forceinline :: dirac_delta
RECURSIVE PURE FUNCTION dirac_delta(E, omega, zeta) RESULT(delta)
  IMPLICIT NONE
  REAL(REAL64) :: delta
  REAL(REAL64), INTENT(IN) :: E, omega, zeta
  delta = zeta / (PI * ((E - omega)**2 + zeta**2))
END FUNCTION dirac_delta

RECURSIVE SUBROUTINE COMPUTE_H_A_1(k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, H_A1)
  IMPLICIT NONE
  REAL(REAL64), INTENT(IN) :: k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B
  INTEGER(INT32), INTENT(IN) :: n_xi
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  COMPLEX(REAL64), INTENT(OUT) :: H_A1(24, 24)

  if (n_xi < 19) error stop "COMPUTE_H_A_1: xi must have at least 19 elements"
  H_A1 = CMPLX(0.0d0, 0.0d0, KIND=REAL64)

#include "H_A_1.proto.f90"
END SUBROUTINE COMPUTE_H_A_1

RECURSIVE SUBROUTINE COMPUTE_H_A_2(k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, H_A2)
  IMPLICIT NONE
  REAL(REAL64), INTENT(IN) :: k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B
  INTEGER(INT32), INTENT(IN) :: n_xi
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  COMPLEX(REAL64), INTENT(OUT) :: H_A2(24, 24)

  if (n_xi < 17) error stop "COMPUTE_H_A_2: xi must have at least 17 elements"
  H_A2 = CMPLX(0.0d0, 0.0d0, KIND=REAL64)

#include "H_A_2.proto.f90"
END SUBROUTINE COMPUTE_H_A_2

SUBROUTINE COMPUTE_H_E(k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, H_E)
  IMPLICIT NONE
  REAL(REAL64), INTENT(IN) :: k_x, k_y, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B
  INTEGER(INT32), INTENT(IN) :: n_xi
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  COMPLEX(REAL64), INTENT(OUT) :: H_E(24, 24)

  if (n_xi < 72) error stop "COMPUTE_H_E: xi must have at least 72 elements"
  H_E = CMPLX(0.0d0, 0.0d0, KIND=REAL64)

#include "H_E.proto.f90"
END SUBROUTINE COMPUTE_H_E

SUBROUTINE COMPUTE_DOS(n_k_points, n_k_points_refined, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, &
                     & zeta, ir_number, xi, DOS_energies, DOS, n_xi, n_energy_points_dos)
  IMPLICIT NONE

  INTEGER(INT32), INTENT(IN) :: n_k_points, n_k_points_refined
  REAL(REAL64), INTENT(IN) :: t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B
  REAL(REAL64), INTENT(IN) :: zeta
  INTEGER(INT32), INTENT(IN) :: ir_number
  INTEGER(INT32), INTENT(IN) :: n_xi, n_energy_points_dos
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  REAL(REAL64), INTENT(IN) :: DOS_energies(n_energy_points_dos)
  REAL(REAL64), INTENT(OUT) :: DOS(n_energy_points_dos, 24)

  COMPLEX(REAL64) :: Hamiltonian(24, 24)
  REAL(REAL64) :: Energies(24)

  REAL(REAL64) :: k1, k2, kx, ky
  REAL(REAL64) :: DOS_local(n_energy_points_dos, 24)

  INTEGER(INT32) :: i, j, n, k, i_ref, j_ref

  REAL(REAL64) :: E_max
  REAL(REAL64) :: dk1, dk2, dk1_ref, dk2_ref

  WRITE (*, *) "Max num threads = ", omp_get_max_threads()

  DOS = 0.0d0
  Hamiltonian = CMPLX(0.0d0, 0.0d0, KIND=REAL64)

  IF (n_k_points <= 0) ERROR STOP "COMPUTE_DOS: n_k_points must be > 0"
  IF (n_k_points_refined < 0) ERROR STOP "COMPUTE_DOS: n_k_points_refined must be >= 0"
  IF (n_energy_points_dos <= 0) ERROR STOP "COMPUTE_DOS: n_energy_points_dos must be > 0"
  IF (ir_number < 1 .or. ir_number > 3) ERROR STOP "COMPUTE_DOS: ir_number must be 1, 2, or 3"

  E_max = MAXVAL(DOS_energies)
  dk1 = K1_MAX / n_k_points
  dk2 = K2_MAX / n_k_points

  dk1_ref = dk1 / (n_k_points_refined + 1)
  dk2_ref = dk2 / (n_k_points_refined + 1)

  !$omp parallel private(i, j, n, k, i_ref, j_ref, k1, k2, kx, ky, Hamiltonian, Energies, DOS_local)
  DOS_local = 0.0d0 !Initialize DOS_local for each thread!
  !$omp do collapse(2)
  DO i = -n_k_points / 2, n_k_points / 2
    DO j = -n_k_points / 2, n_k_points / 2
      k1 = i * dk1
      k2 = j * dk2

      kx = 2.*PI / (SQRT(3.0d0)) * k1
      ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
      Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)

      CALL COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, ir_number, Hamiltonian)

      CALL DIAGONALIZE_HERMITIAN(Hamiltonian, Energies, 24)

      ! If lowest energy is beyond the range we are calculating the DOS for, skip
      ! This might be changed if magnetic field is to be introduced
      IF (MINVAL(ABS(Energies)) > E_max) CYCLE

      !Update DOS for current thread.
      DO n = 1, n_energy_points_dos
        DO k = 1, 24
          DOS_local(n, k) = DOS_local(n, k) + dirac_delta(Energies(k), DOS_energies(n), zeta)
        END DO
      END DO

      !Add a grid refinement here, if lowest energy is in [E_DOS_min, E_DOS_max]
      !Center of the cell
      IF (i .lt. n_k_points / 2 .AND. j .lt. n_k_points / 2) THEN
        DO i_ref = 1, n_k_points_refined
          DO j_ref = 1, n_k_points_refined
            k1 = i * dk1 + i_ref * dk1_ref
            k2 = j * dk2 + j_ref * dk2_ref
            kx = 2.*PI / (SQRT(3.0d0)) * k1
            ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
            Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)

            CALL COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, ir_number, Hamiltonian)

            CALL DIAGONALIZE_HERMITIAN(Hamiltonian, Energies, 24)

            ! If lowest energy is beyond the range we are calculating the DOS for, skip
            ! This might be changed if magnetic field is to be introduced
            IF (MINVAL(ABS(Energies)) > E_max) CYCLE

            !Update DOS for current thread.
            DO n = 1, n_energy_points_dos
              DO k = 1, 24
                DOS_local(n, k) = DOS_local(n, k) + dirac_delta(Energies(k), DOS_energies(n), zeta)
              END DO
            END DO

          END DO
        END DO
      END IF

      IF ((i .lt. n_k_points / 2 .AND. j .lt. n_k_points / 2) .OR. (i .eq. n_k_points / 2)) THEN
        !Left edge without corner
        DO j_ref = 1, n_k_points_refined
          k1 = i * dk1
          k2 = j * dk2 + j_ref * dk2_ref
          kx = 2.*PI / (SQRT(3.0d0)) * k1
          ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
          Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)

          CALL COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, ir_number, Hamiltonian)

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian, Energies, 24)

          ! If lowest energy is beyond the range we are calculating the DOS for, skip
          ! This might be changed if magnetic field is to be introduced
          IF (MINVAL(ABS(Energies)) > E_max) CYCLE

          !Update DOS for current thread.
          DO n = 1, n_energy_points_dos
            DO k = 1, 24
              DOS_local(n, k) = DOS_local(n, k) + dirac_delta(Energies(k), DOS_energies(n), zeta)
            END DO
          END DO

        END DO
      END IF

      IF ((i .lt. n_k_points / 2 .AND. j .lt. n_k_points / 2) .OR. (j .eq. n_k_points / 2)) THEN
        !Bottom edge without corner
        DO i_ref = 1, n_k_points_refined
          k1 = i * dk1 + i_ref * dk1_ref
          k2 = j * dk2
          kx = 2.*PI / (SQRT(3.0d0)) * k1
          ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
          Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
          CALL COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, ir_number, Hamiltonian)

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian, Energies, 24)

          ! If lowest energy is beyond the range we are calculating the DOS for, skip
          ! This might be changed if magnetic field is to be introduced
          IF (MINVAL(ABS(Energies)) > E_max) CYCLE

          !Update DOS for current thread.
          DO n = 1, n_energy_points_dos
            DO k = 1, 24
              DOS_local(n, k) = DOS_local(n, k) + dirac_delta(Energies(k), DOS_energies(n), zeta)
            END DO
          END DO

        END DO
      END IF

    END DO
  END DO
  !$omp end do

  !$omp critical (accumulate_dos)
  DO n = 1, n_energy_points_dos
    DOS(n, :) = DOS(n, :) + DOS_local(n, :)
  END DO
  !$omp end critical (accumulate_dos)

  !$omp end parallel

  ! Normalize the DOS
  IF (MAXVAL(SUM(DOS, DIM=2)) .NE. 0.0d0) DOS = DOS / MAXVAL(SUM(DOS, DIM=2))

END SUBROUTINE COMPUTE_DOS

SUBROUTINE COMPUTE_FREE_ENERGY(n_k_points, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, T,&
                             & ir_number, xi, n_xi, free_energy)
  IMPLICIT NONE

  INTEGER(INT32), INTENT(IN) :: n_k_points
  REAL(REAL64), INTENT(IN) :: t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, T
  INTEGER(INT32), INTENT(IN) :: ir_number
  INTEGER(INT32), INTENT(IN) :: n_xi
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  REAL(REAL64), INTENT(OUT) :: free_energy
  INTEGER(INT32) :: n_triangle, j_phi, i_r
  REAL(REAL64) :: phi_k, r_max, dr, r_k, dphi_k
  REAL(REAL64) :: kx, ky

  COMPLEX(REAL64) :: Hamiltonian(24, 24)
  REAL(REAL64) :: Energies(24)
  INTEGER(INT32) :: i_state, i_xi

  REAL(REAL64) :: Partition_log(24)

  REAL(REAL64) :: free_energy_local

  free_energy = 0.
  dphi_k = (PI / 3.0d0) / n_k_points

  !$omp parallel private(phi_k, r_k, r_max, dr, kx, ky, Hamiltonian, Energies, Partition_log, free_energy_local)
  free_energy_local = 0.
  !$omp do collapse(3)
  DO n_triangle = -N_BZ_SECTIONS / 2, N_BZ_SECTIONS / 2 - 1
    DO j_phi = 0, n_k_points - 1
      DO i_r = 0, n_k_points
        phi_k = n_triangle * (PI / 3.0d0) + j_phi * dphi_k
        r_max = r_max_phi(MOD(ABS(phi_k), PI / 3))
        dr = r_max / n_k_points
        r_k = i_r * dr

        !Transform from graphene reciprocal lattice to kx and ky
        kx = r_k * COS(phi_k)
        ky = r_k * SIN(phi_k)

        CALL COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, &
                                  & xi, n_xi, ir_number, Hamiltonian)

        CALL DIAGONALIZE_HERMITIAN(Hamiltonian, Energies, 24)

        DO i_state = 1, 12
          Partition_log(i_state) = LOG(1.0 + EXP(Energies(i_state) / (k_B * T))) ! Hole states
          Partition_log(12 + i_state) = LOG(1.0 + EXP(-Energies(12 + i_state) / (k_B * T))) ! Electron states
        END DO

        free_energy_local = 0
        DO i_state = 1, 24
          free_energy_local = free_energy_local - (k_B * T) * Partition_log(i_state) * r_k * dr
        END DO

        !Add diagonal terms

        DO i_xi = 1, n_xi
          free_energy_local = free_energy_local + 0.5 * ABS(xi(i_xi)) * r_k * dr! This shall account for condensation energy
        END DO

        !Calculate chemical potential times number of particles

        !$omp critical (free_energy_update)
        free_energy = free_energy + free_energy_local
        !$omp end critical (free_energy_update)
      END DO
    END DO
  END DO
  !$omp end do
  !$omp end parallel

  free_energy = free_energy * JACOBIAN * dphi_k

END SUBROUTINE COMPUTE_FREE_ENERGY

!dir$ attributes forceinline :: COMPUTE_IR_HAMILTONIAN
SUBROUTINE COMPUTE_IR_HAMILTONIAN(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, ir_number, Hamiltonian)
  IMPLICIT NONE
  REAL(REAL64), INTENT(IN) :: kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B
  INTEGER(INT32), INTENT(IN) :: n_xi, ir_number
  COMPLEX(REAL64), INTENT(IN) :: xi(n_xi)
  COMPLEX(REAL64), INTENT(OUT) :: Hamiltonian(24, 24)

  Hamiltonian = CMPLX(0., 0., KIND=REAL64)

  SELECT CASE (ir_number)
  CASE (1)
    CALL COMPUTE_H_A_1(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, Hamiltonian)
  CASE (2)
    CALL COMPUTE_H_A_2(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, Hamiltonian)
  CASE (3)
    CALL COMPUTE_H_E(kx, ky, t_sigma, t_pi, lambda_soc, t_rashba, delta_tri, mu, g, mu_B, B, theta_B, phi_B, xi, n_xi, Hamiltonian)
  CASE DEFAULT
    ERROR STOP "COMPUTE_IR_HAMILTONIAN: ir_number must be 1, 2, or 3"
  END SELECT
END SUBROUTINE COMPUTE_IR_HAMILTONIAN

!dir$ attributes forceinline :: r_max_phi
RECURSIVE FUNCTION r_max_phi(phi) RESULT(r_max)
  IMPLICIT NONE
  REAL(REAL64) :: r_max
  REAL(REAL64), INTENT(IN) :: phi
  r_max = R_K_MAX * SQRT(3.0d0) / (2.0d0 * COS(phi - PI / 6.0d0))
END FUNCTION r_max_phi

END MODULE hamiltonians_sympy
