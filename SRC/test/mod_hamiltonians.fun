test_suite mod_hamiltonians

REAL*8, PARAMETER :: gFactor = 3.0d0
REAL*8, PARAMETER :: muB = 0.5
REAL*8, PARAMETER :: s = 0.5

setup
USE mod_parameters
USE mod_utilities
USE mod_reader
SUBLATTICES = 2
CALL SET_HAMILTONIAN_PARAMS()
end setup

teardown
end teardown

test test_compute_tba_term

REAL*8 :: kx, ky, dk
COMPLEX*16 :: Hamiltonian(DIM, DIM)
INTEGER*4 :: i, j, ik, jk

!Set variables that would have been set by reding input
t_D = .5
t_I = .04
dk = 0.5

DO ik = -2, 2
  DO jk = -2, 2
    Hamiltonian(:, :) = DCMPLX(0., 0.)
    kx = ik * dk
    ky = jk * dk

    CALL COMPUTE_TBA_TERM(Hamiltonian, kx, ky)

    DO i = 1, DIM
      DO j = 1, DIM
        IF (i == 1 .AND. j == 4) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_yz(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_yz(kx, ky)))
        ELSE IF (i == 2 .AND. j == 5) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_zx(kx, ky)))
        ELSE IF (i == 3 .AND. j == 6) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_xy(kx, ky)))
        ELSE IF (i == TBA_DIM + 1 .AND. j == TBA_DIM + 4) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_yz(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_yz(kx, ky)))
        ELSE IF (i == TBA_DIM + 2 .AND. j == TBA_DIM + 5) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_zx(kx, ky)))
        ELSE IF (i == TBA_DIM + 3 .AND. j == TBA_DIM + 6) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(epsilon_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(epsilon_xy(kx, ky)))
          !Nambu space
        ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 4) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_yz(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky))))
        ELSE IF (i == DIM_POSITIVE_K + 2 .AND. j == DIM_POSITIVE_K + 5) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky))))
        ELSE IF (i == DIM_POSITIVE_K + 3 .AND. j == DIM_POSITIVE_K + 6) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky))))
        ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 4) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_yz(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky))))
        ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + 5) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky))))
        ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 3 .AND. j == DIM_POSITIVE_K + TBA_DIM + 6) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(-CONJG(epsilon_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky))))
        ELSE
          assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        END IF
      END DO
    END DO !Check equality

  END DO
END DO !Loop over k

end test

test test_compute_atomic_soc_terms

COMPLEX*16 :: Hamiltonian(DIM, DIM)
INTEGER*4 :: i, j

lambda_SOC = 1.

Hamiltonian(:, :) = DCMPLX(0, 0)
CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian)

DO i = 1, DIM
  DO j = 1, DIM
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

COMPLEX*16 :: Hamiltonian(DIM, DIM)
INTEGER*4 :: i, j

delta_trigonal = 1.
Hamiltonian(:, :) = 0.0d0

CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian)

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

end test

test test_compute_rashba_hopping

COMPLEX*16 :: Hamiltonian(DIM, DIM)
REAL*8 :: kx, ky, dk
INTEGER*4 :: i, j, ik, jk

t_Rashba = 1.
dk = 0.5

DO ik = -2, 2
  DO jk = -2, 2
    kx = ik * dk
    ky = jk * dk

    Hamiltonian(:, :) = 0.0d0
    CALL COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky)

    DO i = 1, DIM
      DO j = 1, DIM
        ! Assert real and imaginary parts of the Hamiltonian matrix
        IF (i == 1 .AND. j == ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_zx(kx, ky)))
        ELSE IF (i == 1 .AND. j == ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_xy(kx, ky)))
        ELSE IF (i == 2 .AND. j == ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_zx(kx, ky)))
        ELSE IF (i == 2 .AND. j == ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_zx_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_zx_xy(kx, ky)))
        ELSE IF (i == 3 .AND. j == ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_xy(kx, ky)))
        ELSE IF (i == 3 .AND. j == ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_zx_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_zx_xy(kx, ky)))
        ELSE IF (i == 1 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_zx(kx, ky)))
        ELSE IF (i == 1 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_yz_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_yz_xy(kx, ky)))
        ELSE IF (i == 2 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_zx(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_zx(kx, ky)))
        ELSE IF (i == 2 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(rashba_zx_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(rashba_zx_xy(kx, ky)))
        ELSE IF (i == 3 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_yz_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_yz_xy(kx, ky)))
        ELSE IF (i == 3 + TBA_DIM .AND. j == TBA_DIM + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(rashba_zx_xy(kx, ky)))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(rashba_zx_xy(kx, ky)))
        ELSE IF (i == 1 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_zx(-kx, -ky))))
        ELSE IF (i == 1 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_xy(-kx, -ky))))
        ELSE IF (i == 2 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_zx(-kx, -ky))))
        ELSE IF (i == 2 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_zx_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_zx_xy(-kx, -ky))))
        ELSE IF (i == 3 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_xy(-kx, -ky))))
        ELSE IF (i == 3 + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_zx_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_zx_xy(-kx, -ky))))
        ELSE IF (i == 1 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_zx(-kx, -ky))))
        ELSE IF (i == 1 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_yz_xy(-kx, -ky))))
        ELSE IF (i == 2 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_zx(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_zx(-kx, -ky))))
        ELSE IF (i == 2 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_zx_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), -AIMAG(CONJG(rashba_zx_xy(-kx, -ky))))
        ELSE IF (i == 3 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_yz_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_yz_xy(-kx, -ky))))
        ELSE IF (i == 3 + TBA_DIM + DIM_POSITIVE_K .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(CONJG(rashba_zx_xy(-kx, -ky))))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(CONJG(rashba_zx_xy(-kx, -ky))))
        ELSE
          assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
        END IF

      END DO
    END DO
  END DO
END DO

CALL COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky)

end test

test test_compute_layer_potential
COMPLEX*16 :: Hamiltonian(DIM, DIM)
INTEGER*4 :: i, j
ALLOCATE (V_layer(SUBLATTICES))

Hamiltonian = 0.0d0
V_layer = (/1.0, 3.0/)

CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian)

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

DEALLOCATE (V_layer)
end test

test test_compute_zeeman
COMPLEX*16 :: Hamiltonian(DIM, DIM)
INTEGER*4 :: i, j, nambu_sign, spin_sign
REAL*8 :: B(3) = [1, 2, 3]

Hamiltonian = 0.0d0
CALL COMPUTE_ZEEMAN(B, Hamiltonian)

DO i = 1, DIM
  nambu_sign = (-1)**(INT(i - 1) / DIM_POSITIVE_K)
  spin_sign = (-1)**(INT(i - 1) / TBA_DIM)
  DO j = 1, DIM
    IF (i == j) THEN
      !B_z terms
      assert_real_equal(REAL(Hamiltonian(i, j)), gFactor * muB * s * nambu_sign * spin_sign * B(3))
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    ELSE IF ( ((i .LE. TBA_DIM) .OR. ((i .GE. DIM_POSITIVE_K + 1) .AND. (i .LE. TBA_DIM + DIM_POSITIVE_K)))&
               .AND. j == i + TBA_DIM) THEN
      !B_x term
      assert_real_equal(REAL(Hamiltonian(i, j)), gFactor * muB * s * nambu_sign * spin_sign * B(1))
      !B_y term
      assert_real_equal(AIMAG(Hamiltonian(i, j)), -gFactor * muB * s * nambu_sign * spin_sign * B(2))
    ELSE
      assert_real_equal(REAL(Hamiltonian(i, j)), 0.0d0)
      assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
    END IF
  END DO
END DO
end test

test test_compute_orbital_magnetic_coupling
COMPLEX*16 :: Hamiltonian(DIM, DIM)
REAL*8, PARAMETER :: muB = 0.5
INTEGER*4 :: i, j, nambu_sign, spin_sign
REAL*8 :: B(3) = [1, 9, 18]
COMPLEX*16 :: elem
COMPLEX*16, PARAMETER :: c_zero = (0.0d0, 0.0d0) ! So that compiler does not complain about type mismatch between 0 and imag
COMPLEX*16, PARAMETER :: L_x(3,3) = TRANSPOSE(RESHAPE([c_zero, c_zero, imag, &
                                                        c_zero, c_zero, imag, &
                                                        -imag, -imag, c_zero], &
                                                        [3,3])) / SQRT(2.0d0)
COMPLEX*16, PARAMETER :: L_y(3,3) = TRANSPOSE(RESHAPE([c_zero, -2*imag, -imag, &
                                                        2*imag, c_zero, imag, &
                                                        imag, -imag, c_zero], &
                                                        [3,3])) / SQRT(6.0d0)
COMPLEX*16, PARAMETER :: L_z(3,3)= TRANSPOSE(RESHAPE([c_zero, -imag, imag, &
                                                        imag, c_zero, -imag, &
                                                        -imag, imag, c_zero], &
                                                        [3,3])) / SQRT(3.0d0)

Hamiltonian = 0.0d0
CALL COMPUTE_ORBITAL_MAGNETIC_COUPLING(B, Hamiltonian)

OPEN(unit=9, file="../../OutputData/OrbitalMagnetic.dat", FORM="FORMATTED", ACTION="WRITE")
DO i = 1, DIM
  WRITE(9,'(24F8.2)') REAL(Hamiltonian(i, :))
  WRITE(9,*)
END DO
WRITE(9,*)
WRITE(9,*)
DO i = 1, DIM
  WRITE(9,'(24F8.2)') AIMAG(Hamiltonian(i, :))
  WRITE(9,*)
END DO
CLOSE(9)

DO i = 1, DIM
  nambu_sign = (-1)**(INT(i - 1) / DIM_POSITIVE_K)
  DO j = 1, DIM
    IF ((MOD(i, ORBITALS) .NE. 0)) THEN
      IF (MOD(i, ORBITALS) .EQ. 1) THEN
        IF (j .EQ. i + 1) THEN
          elem = muB * B(1) * L_x(1,2) * nambu_sign !B_x coupling
          elem = elem + muB * B(2) * L_y(1,2) * nambu_sign !B_y coupling
          elem = elem + muB * B(3) * L_z(1,2) * nambu_sign !B_z coupling
          !Check
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(elem))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(elem))
        ELSE IF (j .EQ. i + 2) THEN
          elem = muB * B(1) * L_x(1,3) * nambu_sign !B_x coupling
          elem = elem + muB * B(2) * L_y(1,3) * nambu_sign !B_y coupling
          elem = elem + muB * B(3) * L_z(1,3) * nambu_sign !B_z coupling
          assert_real_equal(REAL(Hamiltonian(i, j)), REAL(elem))
          assert_real_equal(AIMAG(Hamiltonian(i, j)), AIMAG(elem))
        END IF
      ELSE IF (MOD(i, ORBITALS .EQ. 2)) THEN
        IF (j .EQ. i + 2) THEN
          elem = muB * B(1) * L_x(2,3) * nambu_sign !B_x coupling
          elem = elem + muB * B(2) * L_y(2,3) * nambu_sign !B_y coupling
          elem = elem + muB * B(3) * L_z(2,3) * nambu_sign !B_z coupling
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
end test


end test_suite
