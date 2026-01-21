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

MODULE utilities
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
USE types
IMPLICIT NONE
CONTAINS

SUBROUTINE DIAGONALIZE_HERMITIAN(Hamiltonian, Eigenvalues, N)
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

SUBROUTINE DIAGONALIZE_GENERALIZED(Hamiltonian, Eigenvalues, U_transformation, N)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: N
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(N, N)
  COMPLEX(REAL64), INTENT(INOUT) :: U_transformation(N, N)
  REAL(REAL64), INTENT(OUT) :: Eigenvalues(N)
  COMPLEX(REAL64), ALLOCATABLE :: W(:)
  COMPLEX(REAL64), ALLOCATABLE :: VL(:, :)
  COMPLEX(REAL64), ALLOCATABLE :: WORK(:)
  REAL(REAL64), ALLOCATABLE :: RWORK(:)
  INTEGER(INT32) :: LWORK
  INTEGER(INT32) :: INFO, i

  LWORK = 10 * N
  INFO = 0

  ALLOCATE (W(N))
  ALLOCATE (VL(N, N))
  ALLOCATE (WORK(LWORK))
  ALLOCATE (RWORK(2 * N))

  W(:) = 0.
  VL(:, :) = CMPLX(0., 0., KIND=REAL64)
  WORK(:) = 0.
  RWORK(:) = 0.

  CALL ZGEEV('N', 'V', N, Hamiltonian, N, W, VL, N, U_transformation, N,&
  & WORK, LWORK, RWORK, INFO)

  IF (INFO .ne. 0) THEN
    PRINT *, 'ZGEEV INFO ', INFO
    STOP
  END IF

  !Removing phase ambiguity
  ! DO i = 1, N
  !   U_transformation(:, i) = U_transformation(:, i) * CONJG(U_transformation(1, i)) / ABS(U_transformation(1, i))
  ! END DO
  Eigenvalues(:) = REAL(W(:))

  DEALLOCATE (W)
  DEALLOCATE (VL)
  DEALLOCATE (WORK)
  DEALLOCATE (RWORK)

END SUBROUTINE DIAGONALIZE_GENERALIZED

PURE RECURSIVE SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, N)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: N
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(N, N)
  INTEGER(INT32) :: i, j
  DO i = 1, N
    DO j = i + 1, N
      Hamiltonian(j, i) = CONJG(Hamiltonian(i, j))
    END DO
  END DO
END SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS

!---------------------------------------------------------------------------------------
!------------------------------ KINETIC TERMS ------------------------------------------
!---------------------------------------------------------------------------------------
!dir$ attributes forceinline :: epsilon_yz
PURE FUNCTION epsilon_yz(kx, ky, t_D, t_I) RESULT(epsilon)
  COMPLEX(REAL64) :: epsilon
  REAL(REAL64), INTENT(IN) :: kx, ky, t_D, t_I
  epsilon = -t_D * (EXP(imag * ky) + EXP(imag * (SQRT(3.) / 2.*kx - 1./2.*ky))) &
            - t_I * EXP(-imag * (SQRT(3.) / 2.*kx + 1./2.*ky))
END FUNCTION epsilon_yz

!dir$ attributes forceinline :: epsilon_zx
PURE FUNCTION epsilon_zx(kx, ky, t_D, t_I) RESULT(epsilon)
  COMPLEX(REAL64) :: epsilon
  REAL(REAL64), INTENT(IN) :: kx, ky, t_D, t_I
  epsilon = -t_D * (EXP(imag * ky) + EXP(-imag * (SQRT(3.) / 2.*kx + 1./2.*ky))) &
            - t_I * EXP(imag * (SQRT(3.) / 2.*kx - 1./2.*ky))
END FUNCTION epsilon_zx

!dir$ attributes forceinline :: epsilon_xy
PURE FUNCTION epsilon_xy(kx, ky, t_D, t_I) RESULT(epsilon)
  COMPLEX(REAL64) :: epsilon
  REAL(REAL64), INTENT(IN) :: kx, ky, t_D, t_I
  epsilon = -2.*t_D * COS(SQRT(3.) / 2.*kx) * EXP(-imag * 1./2.*ky) - t_I * EXP(imag * ky)
END FUNCTION epsilon_xy

!---------------------------------------------------------------------------------------
!------------------------------ RASHBA TERMS -------------------------------------------
!---------------------------------------------------------------------------------------
!dir$ attributes forceinline :: rashba_yz_xz
PURE FUNCTION rashba_yz_zx(kx, ky, t_Rashba) RESULT(rashba)
  COMPLEX(REAL64) :: rashba
  REAL(REAL64), INTENT(IN) :: kx, ky, t_Rashba
  rashba = 2 * imag * t_Rashba * SIN(-SQRT(3.) / 2.*kx) * EXP(-imag * 1./2.*ky)
END FUNCTION rashba_yz_zx

!dir$ attributes forceinline :: rashba_yz_xy
PURE FUNCTION rashba_yz_xy(kx, ky, t_Rashba) RESULT(rashba)
  COMPLEX(REAL64) :: rashba
  REAL(REAL64), INTENT(IN) :: kx, ky, t_Rashba
  rashba = -t_Rashba * EXP(imag * ky) * (1.-EXP(imag * (-SQRT(3.) / 2.*kx - 3./2.*ky)))
END FUNCTION rashba_yz_xy

!dir$ attributes forceinline :: rashba_zx_xy
PURE FUNCTION rashba_zx_xy(kx, ky, t_Rashba) RESULT(rashba)
  COMPLEX(REAL64) :: rashba
  REAL(REAL64), INTENT(IN) :: kx, ky, t_Rashba
  rashba = -t_Rashba * EXP(imag * ky) * (1.-EXP(imag * (SQRT(3.) / 2.*kx - 3./2.*ky)))
END FUNCTION rashba_zx_xy

!dir$ attributes forceinline :: compute_nearest_even_hopping
PURE RECURSIVE SUBROUTINE COMPUTE_NEAREST_EVEN_HOPPING(Hoppings, kx, ky, t_D, t_I, n)
  IMPLICIT NONE

  REAL(REAL64), INTENT(IN) :: kx, ky, t_D, t_I
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Hoppings(n)
  REAL(REAL64) :: kx_sqrt3_2, ky_1_2
  COMPLEX(REAL64) :: exp_iky ! exp(imag * ky)
  COMPLEX(REAL64) :: exp_pm ! exp(imag * (kx_sqrt3_2 - ky_1_2))
  COMPLEX(REAL64) :: exp_mm ! exp(-imag * (kx_sqrt3_2 + ky_1_2))

  kx_sqrt3_2 = 0.5 * SQRT(3.) * kx
  ky_1_2 = 0.5 * ky
  !Utilizing Euler's formula for performance
  exp_iky = CMPLX(COS(ky), SIN(ky), REAL64)
  exp_pm = CMPLX(COS(kx_sqrt3_2 - ky_1_2), SIN(kx_sqrt3_2 - ky_1_2), REAL64)
  exp_mm = CMPLX(COS(kx_sqrt3_2 + ky_1_2), -SIN(kx_sqrt3_2 + ky_1_2), REAL64)

  !Electrons
  Hoppings(1) = -t_D * (exp_iky + exp_pm) - t_I * exp_mm !yz - yz
  Hoppings(2) = -t_D * (exp_iky + exp_mm) - t_I * exp_pm !zx - zx
  Hoppings(3) = -2.*t_D * COS(kx_sqrt3_2) * EXP(-imag * ky_1_2) - t_I * exp_iky !xy - xy

  !Holes
  Hoppings(4) = -Hoppings(1) !yz - yz
  Hoppings(5) = -Hoppings(2) !zx - zx
  Hoppings(6) = -Hoppings(3) !xy - xy
END SUBROUTINE COMPUTE_NEAREST_EVEN_HOPPING

!dir$ attributes forceinline :: compute_nearest_odd_hopping
PURE RECURSIVE SUBROUTINE COMPUTE_NEAREST_ODD_HOPPING(Hoppings, kx, ky, coupling_energy, n)
  REAL(REAL64), INTENT(IN) :: kx, ky, coupling_energy
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Hoppings(n)

  REAL(REAL64) :: kx_sqrt3_2, ky_1_2
  COMPLEX(REAL64) :: exp_iky ! exp(imag * ky)
  COMPLEX(REAL64) :: exp_minus_iky_1_2 ! exp(imag * ky_1_2)
  COMPLEX(REAL64) :: exp_mm ! exp(-imag * (kx_sqrt3_2 + ky_1_2))
  COMPLEX(REAL64) :: exp_pm ! exp(imag * (kx_sqrt3_2 - ky_1_2))

  kx_sqrt3_2 = 0.5 * SQRT(3.) * kx
  ky_1_2 = 0.5 * ky

  exp_iky = CMPLX(COS(ky), SIN(ky), REAL64)
  exp_minus_iky_1_2 = CMPLX(COS(ky_1_2), -SIN(ky_1_2), REAL64)
  exp_mm = CMPLX(COS(kx_sqrt3_2 + ky_1_2), -SIN(kx_sqrt3_2 + ky_1_2), REAL64)
  exp_pm = CMPLX(COS(kx_sqrt3_2 - ky_1_2), SIN(kx_sqrt3_2 - ky_1_2), REAL64)

  !Electrons
  Hoppings(1) = -coupling_energy * (2.*imag * exp_minus_iky_1_2 * SIN(kx_sqrt3_2)) !yz - zx
  Hoppings(2) = -coupling_energy * (exp_iky - exp_mm) !yz - xy
  Hoppings(3) = -coupling_energy * (exp_iky - exp_pm) !zx - xy
  !Holes
  Hoppings(4) = -Hoppings(1) !yz - zx
  Hoppings(5) = -Hoppings(2) !yz - xy
  Hoppings(6) = -Hoppings(3) !zx - xy
END SUBROUTINE COMPUTE_NEAREST_ODD_HOPPING

!dir$ attributes forceinline :: compute_next_pi_odd_hopping
PURE RECURSIVE SUBROUTINE COMPUTE_NEXT_PI_ODD_HOPPING(Hoppings, kx, ky, coupling_energy, n)
  REAL(REAL64), INTENT(IN) :: kx, ky
  COMPLEX(REAL64), INTENT(IN) :: coupling_energy
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Hoppings(n)

  REAL(REAL64) :: k1, k2, k3
  REAL(REAL64) :: sum_sin, sin1, sin2, sin3

  k1 = -0.5 * SQRT(3.) * kx + 1.5 * ky
  k2 = -0.5 * SQRT(3.) * kx - 1.5 * ky
  k3 = SQRT(3.) * kx

  sin1 = SIN(k1)
  sin2 = SIN(k2)
  sin3 = SIN(k3)
  sum_sin = sin1 + sin2 + sin3
  !Electrons
  Hoppings(1) = -coupling_energy * (sum_sin + sin3) !yz - zx
  Hoppings(2) = coupling_energy * (sum_sin + sin2) !yz - xy
  Hoppings(3) = -coupling_energy * (sum_sin + sin1) !zx - xy
  !Holes
  !The resulting minus sign calculated analitacally
  Hoppings(4) = -Hoppings(1) !yz - zx
  Hoppings(5) = -Hoppings(2) !yz - xy
  Hoppings(6) = -Hoppings(3) !zx - xy
END SUBROUTINE COMPUTE_NEXT_PI_ODD_HOPPING

!dir$ attributes forceinline :: compute_next_sigma_odd_hopping
PURE RECURSIVE SUBROUTINE COMPUTE_NEXT_SIGMA_ODD_HOPPING(Hoppings, kx, ky, coupling_energy, n)
  REAL(REAL64), INTENT(IN) :: kx, ky
  COMPLEX(REAL64), INTENT(IN) :: coupling_energy
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Hoppings(n)

  REAL(REAL64) :: k1, k2, k3
  REAL(REAL64) :: sin1, sin2, sin3

  k1 = -0.5 * SQRT(3.) * kx + 1.5 * ky
  k2 = -0.5 * SQRT(3.) * kx - 1.5 * ky
  k3 = SQRT(3.) * kx

  sin1 = SIN(k1)
  sin2 = SIN(k2)
  sin3 = SIN(k3)
  !Electrons
  Hoppings(1) = coupling_energy * (sin1 + sin2) !yz - zx
  Hoppings(2) = -coupling_energy * (sin1 + sin3) !yz - xy
  Hoppings(3) = coupling_energy * (sin2 + sin3) !zx - xy
  !Holes
  !The resulting minus sign calculated analitacally
  Hoppings(4) = -Hoppings(1) !yz - zx
  Hoppings(5) = -Hoppings(2) !yz - xy
  Hoppings(6) = -Hoppings(3) !zx - xy

END SUBROUTINE COMPUTE_NEXT_SIGMA_ODD_HOPPING

!---------------------------------------------------------------------------------------
!------------------------------ SC PAIRING TERMS ---------------------------------------
!---------------------------------------------------------------------------------------
!dir$ attributes forceinline :: compute_nearest_pairings
PURE RECURSIVE SUBROUTINE COMPUTE_NEAREST_PAIRINGS(Pairings, kx, ky, n)
  REAL(REAL64), INTENT(IN) :: kx, ky
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Pairings(n)

  REAL(REAL64) :: kx_sqrt3_2, ky_1_2
  REAL(REAL64) :: c1, c2, c3, s1, s2, s3

  kx_sqrt3_2 = 0.5 * SQRT(3.) * kx
  ky_1_2 = 0.5 * ky

  c1 = COS(ky)
  s1 = SIN(ky)

  c2 = COS(kx_sqrt3_2 + ky_1_2)
  s2 = -SIN(kx_sqrt3_2 + ky_1_2)

  c3 = COS(-kx_sqrt3_2 + ky_1_2)
  s3 = -SIN(-kx_sqrt3_2 + ky_1_2)

  Pairings(1) = CMPLX(c1, s1, REAL64)
  Pairings(2) = CMPLX(c2, s2, REAL64)
  Pairings(3) = CMPLX(c3, s3, REAL64)
  Pairings(4) = CMPLX(c1, -s1, REAL64)
  Pairings(5) = CMPLX(c2, -s2, REAL64)
  Pairings(6) = CMPLX(c3, -s3, REAL64)

END SUBROUTINE COMPUTE_NEAREST_PAIRINGS

!dir$ attributes forceinline :: compute_next_pairings
PURE RECURSIVE SUBROUTINE COMPUTE_NEXT_PAIRINGS(Pairings, kx, ky, n)
  REAL(REAL64), INTENT(IN) :: kx, ky
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(OUT) :: Pairings(n)

  REAL(REAL64) :: kx_sqrt_3, kx_sqrt3_2, ky_3_2
  REAL(REAL64) :: c1, c2, c3, s1, s2, s3
  kx_sqrt_3 = SQRT(3.) * kx
  kx_sqrt3_2 = 0.5 * SQRT(3.) * kx
  ky_3_2 = 1.5 * ky

  c1 = COS(kx_sqrt_3)
  s1 = SIN(kx_sqrt_3)

  c2 = COS(kx_sqrt3_2 + ky_3_2)
  s2 = SIN(kx_sqrt3_2 + ky_3_2)

  c3 = COS(-kx_sqrt3_2 + ky_3_2)
  s3 = SIN(-kx_sqrt3_2 + ky_3_2)

  Pairings(1) = CMPLX(c1, -s1, REAL64)
  Pairings(2) = CMPLX(c2, -s2, REAL64)
  Pairings(3) = CMPLX(c3, -s3, REAL64)
  Pairings(4) = CMPLX(c1, s1, REAL64)
  Pairings(5) = CMPLX(c2, s2, REAL64)
  Pairings(6) = CMPLX(c3, s3, REAL64)

END SUBROUTINE COMPUTE_NEXT_PAIRINGS

!dir$ attributes forceinline :: pairing_1
PURE FUNCTION pairing_1(ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: ky
  pairing = CMPLX(COS(ky), SIN(ky), REAL64)
END FUNCTION pairing_1

!dir$ attributes forceinline :: pairing_2
PURE FUNCTION pairing_2(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  pairing = CMPLX(COS(SQRT(3.) / 2.*kx + ky / 2.), -SIN(SQRT(3.) / 2.*kx + ky / 2.), REAL64)
END FUNCTION pairing_2

!dir$ attributes forceinline :: pairing_3
PURE FUNCTION pairing_3(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  pairing = CMPLX(COS(-SQRT(3.) / 2.*kx + ky / 2.), -SIN(-SQRT(3.) / 2.*kx + ky / 2.), REAL64)
END FUNCTION pairing_3

!dir$ attributes forceinline :: pairing_nnn_1
PURE FUNCTION pairing_nnn_1(kx) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx
  !pairing = EXP(-imag * (SQRT(3.) * kx))
  pairing = CMPLX(COS(SQRT(3.) * kx), -SIN(SQRT(3.) * kx), REAL64)
END FUNCTION pairing_nnn_1

!dir$ attributes forceinline :: pairing_nnn_2
PURE FUNCTION pairing_nnn_2(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  !pairing = EXP(-imag * (SQRT(3.) / 2.*kx + 3./2.*ky))
  pairing = CMPLX(COS(SQRT(3.) / 2.*kx + 1.5 * ky), -SIN(SQRT(3.) / 2.*kx + 1.5 * ky), REAL64)
END FUNCTION pairing_nnn_2

!dir$ attributes forceinline :: pairing_nnn_3
PURE FUNCTION pairing_nnn_3(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  !pairing = EXP(-imag * (-SQRT(3.) / 2.*kx + 3./2.*ky))
  pairing = CMPLX(COS(-SQRT(3.) / 2.*kx + 1.5 * ky), -SIN(-SQRT(3.) / 2.*kx + 1.5 * ky), REAL64)
END FUNCTION pairing_nnn_3

!dir$ attributes forceinline :: pairing_nnn_4
PURE FUNCTION pairing_nnn_4(kx) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx
  !pairing = EXP(-imag * (-SQRT(3.) * kx))
  pairing = CMPLX(COS(SQRT(3.) * kx), SIN(SQRT(3.) * kx), REAL64)
END FUNCTION pairing_nnn_4

!dir$ attributes forceinline :: pairing_nnn_5
PURE FUNCTION pairing_nnn_5(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  !pairing = EXP(-imag * (-SQRT(3.) / 2.*kx - 3./2.*ky))
  pairing = CMPLX(COS(-SQRT(3.) / 2.*kx - 1.5 * ky), -SIN(-SQRT(3.) / 2.*kx - 1.5 * ky), REAL64)
END FUNCTION pairing_nnn_5

!dir$ attributes forceinline :: pairing_nnn_6
PURE FUNCTION pairing_nnn_6(kx, ky) RESULT(pairing)
  COMPLEX(REAL64) :: pairing
  REAL(REAL64), INTENT(IN) :: kx, ky
  !pairing = EXP(-imag * (SQRT(3.) / 2.*kx - 3./2.*ky))
  pairing = CMPLX(COS(SQRT(3.) / 2.*kx - 1.5 * ky), -SIN(SQRT(3.) / 2.*kx - 1.5 * ky), REAL64)
END FUNCTION pairing_nnn_6

!dir$ attributes forceinline :: fd_distribution
RECURSIVE PURE FUNCTION fd_distribution(E, E_Fermi, T) RESULT(fd)
  IMPLICIT NONE
  REAL(REAL64) :: fd
  REAL(REAL64), INTENT(IN) :: E, E_Fermi, T
  fd = 0.0
  IF (T .ne. 0.0d0) THEN
    fd = 1./(EXP((E - E_Fermi) / (k_B * T)) + 1.)
  ELSE
    IF (E > E_Fermi) THEN
      fd = 0.
    ELSE IF (E == E_Fermi) THEN
      fd = 0.5
    ELSE IF (E < E_Fermi) THEN
      fd = 1.
    END IF
  END IF
END FUNCTION fd_distribution

!dir$ attributes forceinline :: fd_distribution_derivative_E
RECURSIVE PURE FUNCTION fd_distribution_derivative_E(E, E_Fermi, T) RESULT(der)
  IMPLICIT NONE
  REAL(REAL64) :: der
  REAL(REAL64), INTENT(IN) :: E, E_Fermi, T
  REAL(REAL64) :: x
  x = EXP((E - E_Fermi) / (k_B * T))
  der = -x / ((k_B * T) * (x + 1.0d0)**2)
END FUNCTION fd_distribution_derivative_E

!dir$ attributes forceinline :: dirac_delta
RECURSIVE PURE FUNCTION dirac_delta(E, omega, zeta) RESULT(delta)
  IMPLICIT NONE
  REAL(REAL64) :: delta
  REAL(REAL64), INTENT(IN) :: E, omega, zeta
  delta = zeta / (PI * ((E - omega)**2 + zeta**2))
END FUNCTION dirac_delta

!dir$ attributes forceinline :: r_max_phi
RECURSIVE FUNCTION r_max_phi(phi) RESULT(r_max)
  IMPLICIT NONE
  REAL(REAL64) :: r_max
  REAL(REAL64), INTENT(IN) :: phi
  r_max = R_K_MAX * SQRT(3.0d0) / (2.0d0 * COS(phi - PI / 6.0d0))
END FUNCTION r_max_phi

!---------------------------------------------------------------------------------------
!------------------------------ MATHS AND ALGORITHMS -----------------------------------
!---------------------------------------------------------------------------------------
RECURSIVE FUNCTION partition(X, N, low_bound, high_bound) RESULT(pivot_idx)
  IMPLICIT NONE
  INTEGER(INT32) :: pivot_idx
  INTEGER(INT32), INTENT(IN) :: N
  REAL(REAL64), INTENT(INOUT) :: X(N)
  INTEGER(INT32), INTENT(IN) :: low_bound, high_bound

  INTEGER(INT32) :: i, j
  REAL(REAL64) :: pivot, temp

  pivot = X((low_bound + high_bound) / 2)
  i = low_bound - 1
  j = high_bound + 1

  DO
    DO
      i = i + 1
      IF (X(i) >= pivot) EXIT
    END DO

    DO
      j = j - 1
      IF (X(j) <= pivot) EXIT
    END DO

    IF (i >= j) THEN
      pivot_idx = j
      RETURN
    END IF

    temp = X(i)
    X(i) = X(j)
    X(j) = temp
  END DO
END FUNCTION partition

RECURSIVE SUBROUTINE QSORT(X, N, low_bound, high_bound)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: N
  REAL(REAL64), INTENT(INOUT) :: X(N)
  INTEGER(INT32), INTENT(IN) :: low_bound, high_bound

  INTEGER(INT32) :: pivot_index

  IF (low_bound < high_bound) THEN
    pivot_index = partition(X, N, low_bound, high_bound)
    CALL QSORT(X, N, low_bound, pivot_index)
    CALL QSORT(X, N, pivot_index + 1, high_bound)
  END IF
END SUBROUTINE QSORT

FUNCTION kronecker_product(A, B) result(K)
  !! Given two matrices A and B it computes their Kronecker product K = A ⊗ B
  IMPLICIT NONE
  COMPLEX(REAL64), intent(in) :: A(:, :), B(:, :)
  COMPLEX(REAL64) :: K(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2))
  INTEGER(INT32) :: i, j

  do i = 1, size(A, 1)
    do j = 1, size(A, 2)
      K((i - 1) * size(B, 1) + 1:i * size(B, 1), (j - 1) * size(B, 2) + 1:j * size(B, 2)) = A(i, j) * B
    end do
  end do
END FUNCTION kronecker_product

RECURSIVE FUNCTION det(matrix, n) RESULT(determinant)
  !! Given a square matrix it computes its determinant
  IMPLICIT NONE
  COMPLEX(REAL64) :: determinant
  INTEGER(INT32), INTENT(IN) :: n
  COMPLEX(REAL64), INTENT(IN) :: matrix(n, n)
  INTEGER(INT32) :: IPIV(n)
  INTEGER(INT32) :: info, i

  IPIV(:) = 0.0d0

  CALL ZGETRF(n, n, matrix, n, IPIV, info)
  CALL ZLAPMT(.TRUE., n, n, matrix, n, IPIV)
  determinant = 1.0d0
  DO i = 1, n
    determinant = determinant * matrix(i, i)
  END DO
  RETURN

END FUNCTION det

!---------------------- SPARSE MATRIX STORAGE ---------------------------------------------

PURE RECURSIVE SUBROUTINE MATRICIZE_INTERACTION_TENSOR(Tensor, n_index_elems, Matricized_tensor)
  !! Given a tensor V_{ijkl} it returns a 2D matrix representation V_{(ij)(kl)}
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: n_index_elems !! Number of elements a single index can take
  REAL(REAL64), INTENT(IN) :: Tensor(n_index_elems, n_index_elems, n_index_elems, n_index_elems) !! Original tensor
  REAL(REAL64), INTENT(OUT) :: Matricized_tensor(n_index_elems * n_index_elems, n_index_elems * n_index_elems) !! Matricized tensor

  INTEGER(INT32) :: i, j, k, l
  DO i = 1, n_index_elems
    DO j = 1, n_index_elems
      DO k = 1, n_index_elems
        DO l = 1, n_index_elems
          Matricized_tensor(get_matricized_index(i, j, n_index_elems), get_matricized_index(k, l, n_index_elems)) = Tensor(i, j, k, l)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE MATRICIZE_INTERACTION_TENSOR

!TODO: Utilize that the matrix has to be hermitian/symmetric with respect to certain transformations
SUBROUTINE SAVE_SPARSE_MATRIX_IN_CRS(Matrix, Values, Column_indeces, Row_pointers, n_nonzero, matrix_size)
  !! Given a sparse matrix it saves it in a compressed row storage format.
  !! An assumption is that the matrix has to be hermitian.
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: n_nonzero
  INTEGER(INT32), INTENT(IN) :: matrix_size
  REAL(REAL64), INTENT(IN) :: Matrix(matrix_size, matrix_size)
  REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: Values(:)
  INTEGER(INT32), ALLOCATABLE, INTENT(OUT) :: Column_indeces(:)
  INTEGER(INT32), ALLOCATABLE, INTENT(OUT) :: Row_pointers(:)

  INTEGER(INT32) :: i, j, n_current_nonzero_value

  ALLOCATE (Values(n_nonzero))
  ALLOCATE (Column_indeces(n_nonzero))
  ALLOCATE (Row_pointers(matrix_size + 1))

  Row_pointers(1) = 1
  n_current_nonzero_value = 1
  DO i = 1, matrix_size
    DO j = 1, matrix_size
      IF (Matrix(i, j) .NE. CMPLX(0.0d0, 0.0d0)) THEN
        Values(n_current_nonzero_value) = Matrix(i, j)
        Column_indeces(n_current_nonzero_value) = j
        n_current_nonzero_value = n_current_nonzero_value + 1
      END IF
    END DO
    Row_pointers(i + 1) = n_current_nonzero_value
  END DO

  IF (n_current_nonzero_value - 1 .NE. n_nonzero) STOP "Error in saving sparse matrix in CRS format - number of elements nonzero is not correct"
END SUBROUTINE SAVE_SPARSE_MATRIX_IN_CRS

PURE RECURSIVE FUNCTION calculate_number_of_nonzero_elements(Matrix) RESULT(n)
  IMPLICIT NONE
  INTEGER(INT32) :: n
  REAL(REAL64), INTENT(IN) :: Matrix(:, :)

  INTEGER(INT32) :: i, j

  n = 0
  DO i = 1, SIZE(Matrix, 1)
    DO j = 1, SIZE(Matrix, 2)
      IF (Matrix(i, j) .NE. CMPLX(0.0d0, 0.0d0)) THEN
        n = n + 1
      END IF
    END DO
  END DO
END FUNCTION calculate_number_of_nonzero_elements

PURE RECURSIVE FUNCTION get_matricized_index(i, j, n_index_elems) RESULT(index)
  !! Given two indeces of a tensor V_{ijkl} it returns a single index in 2D matrix representation
  !! V_{(ij)(kl)}
  IMPLICIT NONE
  INTEGER(INT32) :: index
  INTEGER(INT32), INTENT(IN) :: i, j !! Indeces from which the combined single index has to be calculated
  INTEGER(INT32), INTENT(IN) :: n_index_elems !! Number of values the i or j index can take.

  index = (i - 1) * n_index_elems + j
END FUNCTION get_matricized_index

PURE RECURSIVE FUNCTION get_dematricized_indeces(matricized_index, n_index_elems) RESULT(Indeces)
  !! Given a single index in 2D matrix representation it returns the indeces of a tensor V_{ijkl}
  !! from V_{(ij)(kl)}
  IMPLICIT NONE
  INTEGER(INT32) :: Indeces(2)
  INTEGER(INT32), INTENT(IN) :: matricized_index !! Single index in 2D matrix representation
  INTEGER(INT32), INTENT(IN) :: n_index_elems !! Number of values the i or j index can take.

  Indeces(1) = (matricized_index - 1) / n_index_elems + 1
  Indeces(2) = matricized_index - (Indeces(1) - 1) * n_index_elems
END FUNCTION get_dematricized_indeces

PURE RECURSIVE FUNCTION get_index_from_degress_of_freedom(Degrees_of_freedom, discretization) RESULT(index)
  !! Given set of quantum numbers returns and index that will be ascribed to this is in the Hamiltonian
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization !! Discretization parameters
  INTEGER(INT32) :: index !! Index in the Hamiltonian
  INTEGER(INT32), INTENT(IN) :: Degrees_of_freedom(4) !! [orb, lat, spin, nambu]
                                                      !! Orbital number [1 ... ORBITALS]
                                                      !! Sublattice number [1 ... SUBLATTICES]
                                                      !! Spin number [1 2]
                                                      !! Electron/hole part [1 2]

  INTEGER(INT32) :: orb, lat, spin, nambu

  orb = Degrees_of_freedom(1)
  lat = Degrees_of_freedom(2)
  spin = Degrees_of_freedom(3)
  nambu = Degrees_of_freedom(4)

  index = orb + (lat - 1) * discretization % ORBITALS + (spin - 1) * discretization % derived % TBA_DIM + &
    & (nambu - 1) * discretization % derived % DIM_POSITIVE_K
END FUNCTION get_index_from_degress_of_freedom

PURE RECURSIVE FUNCTION get_degress_of_freedom_from_index(index, discretization) RESULT(Degrees_of_freedom)
  !! Given an index in the Hamiltonian returns the set of quantum numbers
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization !! Discretization parameters
  INTEGER(INT32) :: Degrees_of_freedom(4) !! [orb, lat, spin, nambu]
                                          !! Orbital number [1 ... ORBITALS]
                                          !! Sublattice number [1 ... SUBLATTICES]
                                          !! Spin number [1 2]
                                          !! Electron/hole part [1 2]
  INTEGER(INT32), INTENT(IN) :: index !! Single index in 2D matrix representation

  Degrees_of_freedom(1) = MOD(index - 1, discretization % ORBITALS) + 1 ! orbital
  Degrees_of_freedom(2) = MOD((index - 1) / discretization % ORBITALS, discretization % SUBLATTICES) + 1 ! sublattice
  Degrees_of_freedom(3) = MOD((index - 1) / discretization % derived % TBA_DIM, SPINS) + 1 ! spin
  Degrees_of_freedom(4) = (index - 1) / discretization % derived % DIM_POSITIVE_K + 1 ! nambu
END FUNCTION get_degress_of_freedom_from_index

! Testing functions
SUBROUTINE ROTATE_HAMILTONIAN_60_DEG(Hamiltonian, DIM)
  !! Rotates the hamiltonian by 60 degrees.
  !! Takse into account orbital, wavevector and spin degrees of freedom
  INTEGER(INT32), INTENT(IN) :: DIM
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(DIM, DIM)

  COMPLEX(REAL64) :: U_orb(3, 3)
  COMPLEX(REAL64) :: U_lat(2, 2)
  COMPLEX(REAL64) :: U_spin(2, 2)
  COMPLEX(REAL64) :: Unity(2, 2)
  COMPLEX(REAL64), ALLOCATABLE :: U_orb_lat(:, :)
  COMPLEX(REAL64), ALLOCATABLE :: U_orb_lat_spin(:, :)
  COMPLEX(REAL64) :: U_orb_lat_spin_nambu(DIM, DIM)

  INTEGER(INT32) :: i, j
  COMPLEX(REAL64), PARAMETER :: c_zero = CMPLX(0.0d0, 0.0d0)
  !Currently only for 3x3 matrices
  !Assuming {yz, zx, xy}
  U_orb = TRANSPOSE(RESHAPE([ &
                            0, 0, 1, &
                            1, 0, 0, &
                            0, 1, 0], [3, 3]))

  ! CONJG_TRANSPOSE(eigvecs) is equivalent to V^{-1} for unitary eigvecs
  !Assuming 2 sublattices
  U_lat = TRANSPOSE(RESHAPE([0, 1, &
                           & 1, 0], [2, 2]))

  U_spin = TRANSPOSE(RESHAPE([EXP(-imag * PI / 6.), c_zero, &
                            & c_zero, EXP(imag * PI / 6.)], [2, 2]))

  Unity = TRANSPOSE(RESHAPE([1, 0, &
                           & 0, 1], [2, 2]))

  U_orb_lat = kronecker_product(U_lat, U_orb)
  U_orb_lat_spin = kronecker_product(U_spin, U_orb_lat)
  U_orb_lat_spin_nambu = c_zero
  DO i = 1, DIM / 2
    DO j = 1, DIM / 2
      U_orb_lat_spin_nambu(i, j) = U_orb_lat_spin(i, j)
      U_orb_lat_spin_nambu(DIM / 2 + i, DIM / 2 + j) = CONJG(U_orb_lat_spin(i, j))
    END DO
  END DO

  Hamiltonian = MATMUL(U_orb_lat_spin_nambu, MATMUL(Hamiltonian, TRANSPOSE(CONJG(U_orb_lat_spin_nambu))))

END SUBROUTINE ROTATE_HAMILTONIAN_60_DEG

END MODULE utilities
