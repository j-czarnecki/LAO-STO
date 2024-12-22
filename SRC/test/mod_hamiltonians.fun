test_suite mod_hamiltonians

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
      Hamiltonian(:,:) = DCMPLX(0., 0.)
      kx = ik*dk
      ky = jk*dk

      CALL COMPUTE_TBA_TERM(Hamiltonian, kx, ky)

      DO i = 1, DIM
        DO j = 1, DIM
          IF (i == 1 .AND. j == 4) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_yz(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_yz(kx, ky)))
          ELSE IF (i == 2 .AND. j == 5) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_zx(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_zx(kx, ky)))
          ELSE IF (i == 3 .AND. j == 6) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_xy(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_xy(kx, ky)))
          ELSE IF (i == TBA_DIM + 1 .AND. j == TBA_DIM + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_yz(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_yz(kx, ky)))
          ELSE IF (i == TBA_DIM + 2 .AND. j == TBA_DIM + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_zx(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_zx(kx, ky)))
          ELSE IF (i == TBA_DIM + 3 .AND. j == TBA_DIM + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(epsilon_xy(kx, ky)))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(epsilon_xy(kx, ky)))
          !Nambu space
          ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_yz(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky))))
          ELSE IF (i == DIM_POSITIVE_K + 2 .AND. j == DIM_POSITIVE_K + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_zx(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky))))
          ELSE IF (i == DIM_POSITIVE_K + 3 .AND. j == DIM_POSITIVE_K + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_xy(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 4) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_yz(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_yz(-kx, -ky))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + 5) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_zx(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_zx(-kx, -ky))))
          ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 3 .AND. j == DIM_POSITIVE_K + TBA_DIM + 6) THEN
            assert_real_equal(REAL(Hamiltonian(i,j)), REAL(-CONJG(epsilon_xy(-kx, -ky))))
            assert_real_equal(AIMAG(Hamiltonian(i,j)), AIMAG(-CONJG(epsilon_xy(-kx, -ky))))
          ELSE
            assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
            assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
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

  Hamiltonian(:,:) = DCMPLX(0, 0)
  CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian)

  DO i = 1, DIM
    DO j = 1, DIM
      IF ((i == 1 .AND. j == 2) .OR. (i == 4 .AND. j == 5)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), lambda_SOC/2.)
      ELSE IF ((i == 7 .AND. j == 8) .OR. (i == 10 .AND. j == 11)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), -lambda_SOC/2.)
      ELSE IF ((i == 1 .AND. j == 9) .OR. (i == 4 .AND. j == 12)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), -lambda_SOC/2.)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      ELSE IF ((i == 2 .AND. j == 9) .OR. (i == 5 .AND. j == 12)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), lambda_SOC/2.)
      ELSE IF ((i == 3 .AND. j == 7) .OR. (i == 6 .AND. j == 10)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), lambda_SOC/2.)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      ELSE IF ((i == 3 .AND. j == 8) .OR. (i == 6 .AND. j == 11)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), -lambda_SOC/2.)
      !Nambu space
      ELSE IF ((i == 13 .AND. j == 14) .OR. (i == 16 .AND. j == 17)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), lambda_SOC/2.)
      ELSE IF ((i == 19 .AND. j == 20) .OR. (i == 22 .AND. j == 23)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), -lambda_SOC/2.)
      ELSE IF ((i == 13 .AND. j == 21) .OR. (i == 16 .AND. j == 24)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), lambda_SOC/2.)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      ELSE IF ((i == 14 .AND. j == 21) .OR. (i == 17 .AND. j == 24)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), lambda_SOC/2.)
      ELSE IF ((i == 15 .AND. j == 19) .OR. (i == 18 .AND. j == 22)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), -lambda_SOC/2.)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      ELSE IF ((i == 15 .AND. j == 20) .OR. (i == 18 .AND. j == 23)) THEN
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), -lambda_SOC/2.)
      ELSE
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      END IF
    END DO
  END DO
end test


test test_compute_trigonal_terms

  COMPLEX*16 :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: i, j

  DELTA_TRI = 1.
  Hamiltonian(:,:) = 0.0d0

  CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian)

  DO i = 1, DIM
    DO j = 1, DIM
      IF (i == 1 .AND. j == 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 .AND. j == 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 .AND. j == 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == ORBITALS + 1 .AND. j == ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == ORBITALS + 1 .AND. j == ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == ORBITALS + 2 .AND. j == ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == 1 + TBA_DIM .AND. j == 2 + TBA_DIM) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 + TBA_DIM .AND. j == 3 + TBA_DIM) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 + TBA_DIM .AND. j == 3 + TBA_DIM) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == 1 + TBA_DIM + ORBITALS .AND. j == 2 + TBA_DIM + ORBITALS) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 1 + TBA_DIM + ORBITALS .AND. j == 3 + TBA_DIM + ORBITALS) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == 2 + TBA_DIM + ORBITALS .AND. j == 3 + TBA_DIM + ORBITALS) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)


      ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + 1 .AND. j == DIM_POSITIVE_K + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + 2 .AND. j == DIM_POSITIVE_K + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + ORBITALS + 2 .AND. j == DIM_POSITIVE_K + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)

      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 1 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE IF (i == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 2 .AND. j == DIM_POSITIVE_K + TBA_DIM + ORBITALS + 3) THEN
          assert_real_equal(REAL(Hamiltonian(i, j)), -DELTA_TRI / 2.0)
          assert_real_equal(AIMAG(Hamiltonian(i, j)), 0.0d0)
      ELSE
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
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
      kx = ik*dk
      ky = jk*dk

      Hamiltonian(:,:) = 0.0d0
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
              assert_real_equal(REAL(Hamiltonian(i, j)), -REAL(CONJG(rashba_yz_xy(-kx,-ky))))
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
  ALLOCATE(V_layer(SUBLATTICES))

  Hamiltonian = 0.0d0
  V_layer = (/ 1.0, 3.0 /)

  CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian)

  DO i = 1, DIM
    DO j = 1, DIM
      IF (i == j) THEN
        IF ((i - 1) / ORBITALS == 0 .OR. (i - 1) / ORBITALS == 2 ) THEN
          assert_real_equal(REAL(Hamiltonian(i,j)), V_layer(1))
          assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 1 .OR. (i - 1) / ORBITALS == 3 ) THEN
          assert_real_equal(REAL(Hamiltonian(i,j)), V_layer(2))
          assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 4 .OR. (i - 1) / ORBITALS == 6 ) THEN
          assert_real_equal(REAL(Hamiltonian(i,j)), -V_layer(1))
          assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
        ELSE IF ((i - 1) / ORBITALS == 5 .OR. (i - 1) / ORBITALS == 7 ) THEN
          assert_real_equal(REAL(Hamiltonian(i,j)), -V_layer(2))
          assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
        END IF
      ELSE
        assert_real_equal(REAL(Hamiltonian(i,j)), 0.0d0)
        assert_real_equal(AIMAG(Hamiltonian(i,j)), 0.0d0)
      END IF
    END DO
  END DO

  DEALLOCATE(V_layer)
end test

end test_suite