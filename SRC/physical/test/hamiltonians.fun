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

test_suite hamiltonians

REAL(REAL64), PARAMETER :: muB = 0.5
REAL(REAL64), PARAMETER :: s = 0.5

TYPE(sc_input_params_t) :: sc_input

setup
INTEGER(INT32) :: sublats = 2
INTEGER(INT32) :: n_subbands = 2
CALL SET_HAMILTONIAN_PARAMS(sublats, n_subbands, sc_input % discretization)

end setup

teardown
end teardown

test test_compute_tba_term

REAL(REAL64) :: kx, ky, dk
COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
INTEGER(INT32) :: i, j, ik, jk

!Set variables that would have been set by reding input
REAL(REAL64) :: t_D = .5
REAL(REAL64) :: t_I = .04
dk = 0.5

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)
  DO ik = -2, 2
    DO jk = -2, 2
      Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
      kx = ik * dk
      ky = jk * dk

      CALL COMPUTE_TBA_TERM(Hamiltonian, kx, ky, t_D, t_I, sc_input % discretization)

      DO i = 1, DIM
        DO j = 1, DIM
          IF (i == 1 .AND. j == 4) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_yz(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_yz(kx, ky, t_D, t_I)))
          ELSE IF (i == 2 .AND. j == 5) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_zx(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_zx(kx, ky, t_D, t_I)))
          ELSE IF (i == 3 .AND. j == 6) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_xy(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_xy(kx, ky, t_D, t_I)))
          ELSE IF (i == TBA_DIM + 1 .AND. j == TBA_DIM + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_yz(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_yz(kx, ky, t_D, t_I)))
          ELSE IF (i == TBA_DIM + 2 .AND. j == TBA_DIM + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_zx(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_zx(kx, ky, t_D, t_I)))
          ELSE IF (i == TBA_DIM + 3 .AND. j == TBA_DIM + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_xy(kx, ky, t_D, t_I)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_xy(kx, ky, t_D, t_I)))
            !Nambu space
          ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_yz(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky, t_D, t_I))))
          ELSE IF (i == DIM_POSITIVE_K + 2 .AND. j == DIM_POSITIVE_K + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_zx(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky, t_D, t_I))))
          ELSE IF (i == DIM_POSITIVE_K + 3 .AND. j == DIM_POSITIVE_K + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_xy(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky, t_D, t_I))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_yz(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky, t_D, t_I))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_zx(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky, t_D, t_I))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 3 .AND. j == DIM_POSITIVE_K + TBA_DIM + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_xy(-kx, -ky, t_D, t_I))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky, t_D, t_I))))
          ELSE
            assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
            assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
          END IF
        END DO
      END DO !Check equality

    END DO
  END DO !Loop over k
END ASSOCIATE
end test

test test_compute_atomic_soc_terms

COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
INTEGER(INT32) :: i, j

REAL(REAL64) :: lambda_SOC = 1.

Hamiltonian(:, :) = CMPLX(0, 0, KIND=REAL64)
CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, lambda_SOC, sc_input % discretization)

DO i = 1, sc_input % discretization % derived % DIM
  DO j = 1, sc_input % discretization % derived % DIM
    IF ((i == 1 .AND. j == 2) .OR. (i == 4 .AND. j == 5)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), lambda_SOC / 2.)
    ELSE IF ((i == 7 .AND. j == 8) .OR. (i == 10 .AND. j == 11)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), -lambda_SOC / 2.)
    ELSE IF ((i == 1 .AND. j == 9) .OR. (i == 4 .AND. j == 12)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), -lambda_SOC / 2.)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    ELSE IF ((i == 2 .AND. j == 9) .OR. (i == 5 .AND. j == 12)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), lambda_SOC / 2.)
    ELSE IF ((i == 3 .AND. j == 7) .OR. (i == 6 .AND. j == 10)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), lambda_SOC / 2.)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    ELSE IF ((i == 3 .AND. j == 8) .OR. (i == 6 .AND. j == 11)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), -lambda_SOC / 2.)
      !Nambu space
    ELSE IF ((i == 13 .AND. j == 14) .OR. (i == 16 .AND. j == 17)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), lambda_SOC / 2.)
    ELSE IF ((i == 19 .AND. j == 20) .OR. (i == 22 .AND. j == 23)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), -lambda_SOC / 2.)
    ELSE IF ((i == 13 .AND. j == 21) .OR. (i == 16 .AND. j == 24)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), lambda_SOC / 2.)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    ELSE IF ((i == 14 .AND. j == 21) .OR. (i == 17 .AND. j == 24)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), lambda_SOC / 2.)
    ELSE IF ((i == 15 .AND. j == 19) .OR. (i == 18 .AND. j == 22)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), -lambda_SOC / 2.)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    ELSE IF ((i == 15 .AND. j == 20) .OR. (i == 18 .AND. j == 23)) THEN
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), -lambda_SOC / 2.)
    ELSE
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    END IF
  END DO
END DO
end test

test test_compute_trigonal_terms

COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
INTEGER(INT32) :: i, j

REAL(REAL64) :: delta_trigonal = 1.
Hamiltonian(:, :) = 0.0d0

CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian, delta_trigonal, sc_input % discretization)

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  DO i = 1, DIM
    DO j = 1, DIM
      IF (i == 1 .AND. j == 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 .AND. j == 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 .AND. j == 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == ORBITALS + 1 .AND. j == ORBITALS + 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == ORBITALS + 1 .AND. j == ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == ORBITALS + 2 .AND. j == ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == 1 + TBA_DIM .AND. j == 2 + TBA_DIM) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 + TBA_DIM .AND. j == 3 + TBA_DIM) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 + TBA_DIM .AND. j == 3 + TBA_DIM) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == 1 + TBA_DIM + ORBITALS .AND. j == 2 + TBA_DIM + ORBITALS) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 + TBA_DIM + ORBITALS .AND. j == 3 + TBA_DIM + ORBITALS) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 + TBA_DIM + ORBITALS .AND. j == 3 + TBA_DIM + ORBITALS) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + 2 .AND. j == DIM_POSITIVE_K + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 2 .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
        assert_real_equal(REAL(Hamiltonian(i, j)), -delta_trigonal / 2.0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE
        assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      END IF
    END DO
  END DO
END ASSOCIATE
end test

test test_compute_rashba_hopping
USE writers
COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
REAL(REAL64) :: kx, ky, dk
INTEGER(INT32) :: i, j, ik, jk

REAL(REAL64) :: t_Rashba = 1.
dk = 0.5

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  DO ik = -2, 2
    DO jk = -2, 2
      kx = ik * dk
      ky = jk * dk

      Hamiltonian(:, :) = 0.0d0
      CALL COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky, t_Rashba, sc_input % discretization)
      CALL PRINT_HAMILTONIAN(Hamiltonian, sc_input % discretization % derived % DIM, "H_rashba")

      DO i = 1, DIM
        DO j = 1, DIM
          ! Assert real and imaginary parts of the Hamiltonian matrix
          IF (i == 1 .AND. j == ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_zx(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_zx(kx, ky, t_Rashba)))
          ELSE IF (i == 1 .AND. j == ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 2 .AND. j == ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_zx(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_zx(kx, ky, t_Rashba)))
          ELSE IF (i == 2 .AND. j == ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_zx_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_zx_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 3 .AND. j == ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 3 .AND. j == ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_zx_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_zx_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 1 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_zx(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_zx(kx, ky, t_Rashba)))
          ELSE IF (i == 1 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 2 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_zx(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_zx(kx, ky, t_Rashba)))
          ELSE IF (i == 2 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_zx_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_zx_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 3 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 3 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_zx_xy(kx, ky, t_Rashba)))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_zx_xy(kx, ky, t_Rashba)))
          ELSE IF (i == 1 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
          ELSE IF (i == 1 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 2 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
          ELSE IF (i == 2 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 3 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 3 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 1 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
          ELSE IF (i == 1 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 2 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))))
          ELSE IF (i == 2 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 3 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))))
          ELSE IF (i == 3 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))))
          ELSE
            assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
            assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
          END IF

        END DO
      END DO
    END DO
  END DO

END ASSOCIATE
end test

test test_compute_fermi_energy
!This test has bad name, it is not dedicated to some particular hamiltonian, rather checking the whole symmetry
USE writers
USE utilities
COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
COMPLEX(REAL64) :: Hamiltonian_rotated(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)

REAL(REAL64) :: kx, ky, kx_rotated, ky_rotated
REAL(REAL64) :: r_k = 1
REAL(REAL64) :: phi_k = 2
REAL(REAL64) :: rotation_angle = PI / 3
INTEGER(INT32) :: n_pi_3_rotations = 2
INTEGER(INT32) :: i, j, ik, jk

REAL(REAL64) :: t_Rashba = 1.
REAL(REAL64) :: t_D = .5
REAL(REAL64) :: t_I = .04
REAL(REAL64) :: lambda_soc = 1.
COMPLEX(REAL64) :: trace = 0.

kx = r_k * COS(phi_k)
ky = r_k * SIN(phi_k)

kx_rotated = r_k * COS(phi_k + n_pi_3_rotations * rotation_angle)
ky_rotated = r_k * SIN(phi_k + n_pi_3_rotations * rotation_angle)

Hamiltonian = CMPLX(0., 0., KIND=REAL64)
Hamiltonian_rotated = CMPLX(0., 0., KIND=REAL64)

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  ! Compute initial hamiltonian
  CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, lambda_soc, sc_input % discretization)
  !CALL COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky, t_Rashba, sc_input % discretization)
  !CALL COMPUTE_TBA_TERM(Hamiltonian, kx, ky, t_D, t_I, sc_input % discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, sc_input % discretization % derived % DIM)
  CALL PRINT_HAMILTONIAN(Hamiltonian, sc_input % discretization % derived % DIM, "H_rashba_original")

  ! Compute rotated hamiltonian
  CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_rotated, lambda_soc, sc_input % discretization)
  !CALL COMPUTE_RASHBA_HOPPING(Hamiltonian_rotated, kx_rotated, ky_rotated, t_Rashba, sc_input % discretization)
  !CALL COMPUTE_TBA_TERM(Hamiltonian_rotated, kx_rotated, ky_rotated, t_D, t_I, sc_input % discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_rotated, sc_input % discretization % derived % DIM)
  CALL PRINT_HAMILTONIAN(Hamiltonian_rotated, sc_input % discretization % derived % DIM, "H_rashba_rotated")

  ! Check that when I rotate by pi/3, the Hamiltonian is invariant
  DO i = 1, n_pi_3_rotations
    CALL ROTATE_HAMILTONIAN_60_DEG(Hamiltonian, sc_input % discretization % derived % DIM)
  END DO
  CALL PRINT_HAMILTONIAN(Hamiltonian, sc_input % discretization % derived % DIM, "H_rashba_original_rotated")

END ASSOCIATE
end test

test test_compute_layer_potential
COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
INTEGER(INT32) :: i, j
REAL(REAL64), ALLOCATABLE :: V_layer(:)
ALLOCATE (V_layer(sc_input % discretization % SUBLATTICES))

Hamiltonian = 0.0d0
V_layer = (/1.0, 3.0/)

CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian, V_layer, sc_input % discretization)

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  DO i = 1, DIM
    DO j = 1, DIM
      IF (i == j) THEN
        IF ((i - 1) / ORBITALS == 0 .OR. (i - 1) / ORBITALS == 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), V_layer(1))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 1 .OR. (i - 1) / ORBITALS == 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), V_layer(2))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 4 .OR. (i - 1) / ORBITALS == 6) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -V_layer(1))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 5 .OR. (i - 1) / ORBITALS == 7) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -V_layer(2))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        END IF
      ELSE
        assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      END IF
    END DO
  END DO
END ASSOCIATE
DEALLOCATE (V_layer)
end test

test test_compute_zeeman
COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
INTEGER(INT32) :: i, j, nambu_sign, spin_sign
REAL(REAL64) :: B(3) = [1, 2, 3]
REAL(REAL64) :: g_factor = 6.0

Hamiltonian = 0.0d0
CALL COMPUTE_ZEEMAN(B, g_factor, Hamiltonian, sc_input % discretization)

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  DO i = 1, DIM
    nambu_sign = (-1)**(INT(i - 1) / DIM_POSITIVE_K)
    spin_sign = (-1)**(INT(i - 1) / TBA_DIM)
    DO j = 1, DIM
      IF (i == j) THEN
        !B_z terms
        assert_real_equal(REAL(Hamiltonian(i, j)), g_factor * muB * s * nambu_sign * spin_sign * B(3))
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (((i .LE. TBA_DIM) .OR. ((i .GE. DIM_POSITIVE_K + 1) .AND. (i .LE. TBA_DIM + DIM_POSITIVE_K))) &
               .AND. j == i + TBA_DIM) THEN
        !B_x term
        assert_real_equal(REAL(Hamiltonian(i, j)), g_factor * muB * s * nambu_sign * spin_sign * B(1))
        !B_y term
        assert_real_equal(AIMAG(Hamiltonian(i, j)), -g_factor * muB * s * nambu_sign * spin_sign * B(2))
      ELSE
        assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      END IF
    END DO
  END DO
END ASSOCIATE
end test

test test_compute_orbital_magnetic_coupling
USE writers

COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
REAL(REAL64), PARAMETER :: muB = 0.5
INTEGER(INT32) :: i, j, nambu_sign, spin_sign
REAL(REAL64) :: B(3) = [1, 9, 18]
COMPLEX(REAL64) :: elem
COMPLEX(REAL64), PARAMETER :: c_zero = (0.0d0, 0.0d0) ! So that compiler does not complain about type mismatch between 0 and imag
COMPLEX(REAL64), PARAMETER :: L_x(3, 3) = TRANSPOSE(RESHAPE([c_zero, c_zero, imag, &
                                                             c_zero, c_zero, imag, &
                                                             -imag, -imag, c_zero], &
                                                            [3, 3])) / SQRT(2.0d0)
COMPLEX(REAL64), PARAMETER :: L_y(3, 3) = TRANSPOSE(RESHAPE([c_zero, -2 * imag, -imag, &
                                                             2 * imag, c_zero, imag, &
                                                             imag, -imag, c_zero], &
                                                            [3, 3])) / SQRT(6.0d0)
COMPLEX(REAL64), PARAMETER :: L_z(3, 3) = TRANSPOSE(RESHAPE([c_zero, -imag, imag, &
                                                             imag, c_zero, -imag, &
                                                             -imag, imag, c_zero], &
                                                            [3, 3])) / SQRT(3.0d0)

Hamiltonian = 0.0d0
CALL COMPUTE_ORBITAL_MAGNETIC_COUPLING(B, Hamiltonian, sc_input % discretization)
CALL PRINT_HAMILTONIAN(Hamiltonian, sc_input % discretization % derived % DIM, "H_L_dot_B")

ASSOCIATE (DIM => sc_input % discretization % derived % DIM, &
           DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
           TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
           LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS, &
           ORBITALS => sc_input % discretization % ORBITALS)

  DO i = 1, DIM
    nambu_sign = (-1)**(INT(i - 1) / DIM_POSITIVE_K)
    DO j = 1, DIM
      IF ((MOD(i, ORBITALS) .NE. 0)) THEN
        IF (MOD(i, ORBITALS) .EQ. 1) THEN
          IF (j .EQ. i + 1) THEN
            elem = muB * B(1) * L_x(1, 2) * nambu_sign !B_x coupling
            elem = elem + muB * B(2) * L_y(1, 2) * nambu_sign !B_y coupling
            elem = elem + muB * B(3) * L_z(1, 2) * nambu_sign !B_z coupling
            !Check
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(elem))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(elem))
          ELSE IF (j .EQ. i + 2) THEN
            elem = muB * B(1) * L_x(1, 3) * nambu_sign !B_x coupling
            elem = elem + muB * B(2) * L_y(1, 3) * nambu_sign !B_y coupling
            elem = elem + muB * B(3) * L_z(1, 3) * nambu_sign !B_z coupling
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(elem))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(elem))
          END IF
        ELSE IF (MOD(i, ORBITALS) .EQ. 2) THEN
          IF (j .EQ. i + 1) THEN
            elem = muB * B(1) * L_x(2, 3) * nambu_sign !B_x coupling
            elem = elem + muB * B(2) * L_y(2, 3) * nambu_sign !B_y coupling
            elem = elem + muB * B(3) * L_z(2, 3) * nambu_sign !B_z coupling
            assert_real_equal(REAL(Hamiltonian(i, j)), REAL(elem))
            assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(elem))
          END IF
        END IF
      ELSE
        assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      END IF
    END DO
  END DO
END ASSOCIATE
end test

test test_compute_sc
USE writers

COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM, sc_input % discretization % derived % DIM)
COMPLEX(REAL64) :: Gamma_SC(sc_input % discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, sc_input % discretization % derived % LAYER_COUPLINGS) !! Superconducting energies
REAL(REAL64) :: kx, ky, dk
INTEGER(INT32) :: i, j, ik, jk

REAL(REAL64) :: gamma_nearest = 1.
REAL(REAL64) :: gamma_next = 3.

dk = 0.

Gamma_SC = 0.
Gamma_SC(:, :N_NEIGHBOURS, 1, 1, :) = CMPLX(gamma_nearest, 0., KIND=REAL64)
Gamma_SC(:, :N_NEIGHBOURS, 2, 2, :) = CMPLX(-gamma_nearest, 0., KIND=REAL64)
!coupling for next nearest neighbours
Gamma_SC(:, (N_NEIGHBOURS + 1):, 1, 1, :) = CMPLX(gamma_next, 0., KIND=REAL64)
Gamma_SC(:, (N_NEIGHBOURS + 1):, 2, 2, :) = CMPLX(-gamma_next, 0., KIND=REAL64)

DO ik = -2, 2
  DO jk = -2, 2

    Hamiltonian = 0.0d0
    kx = ik * dk
    ky = jk * dk

    CALL COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, sc_input % discretization)

    CALL PRINT_HAMILTONIAN(Hamiltonian, sc_input % discretization % derived % DIM, "H_SC")
  END DO
END DO

end test

end test_suite
