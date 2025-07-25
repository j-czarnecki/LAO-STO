MODULE mod_hamiltonians
USE mod_utilities
USE mod_parameters
USE mod_reader
IMPLICIT NONE
CONTAINS


RECURSIVE SUBROUTINE COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian)
  !! Computes all terms that do not depend on k, including complex conjugate elements
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !! Hamiltonian of the system that is to be filled
  CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian(:, :))
  CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian(:, :))
  CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian(:, :))
  CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian(:, :))
  CALL COMPUTE_TETRAGONAL_STRAIN(Hamiltonian(:, :))
  CALL COMPUTE_FERMI_ENERGY(Hamiltonian(:, :))
  CALL COMPUTE_ZEEMAN(B_field, Hamiltonian(:, :))
  CALL COMPUTE_ORBITAL_MAGNETIC_COUPLING(B_field, Hamiltonian(:,:))
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle
END SUBROUTINE COMPUTE_K_INDEPENDENT_TERMS

RECURSIVE SUBROUTINE COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky)
  !! Computes all terms that depend on k, excluding complex conjugate elements
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !! Hamiltonian of the system that is to be filled
  REAL*8, INTENT(INOUT) :: kx !! Wavevector in X direction
  REAL*8, INTENT(INOUT) :: ky !! Wavevector in Y direction
  CALL COMPUTE_TBA_TERM(Hamiltonian(:, :), kx, ky)
  CALL COMPUTE_TI1_TI2(Hamiltonian(:, :), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
  CALL COMPUTE_H_PI(Hamiltonian(:, :), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
  CALL COMPUTE_H_SIGMA(Hamiltonian(:, :), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
  CALL COMPUTE_RASHBA_HOPPING(Hamiltonian(:, :), kx, ky) !This is adapted from KTaO_3, see: PRB, 103, 035115
END SUBROUTINE COMPUTE_K_DEPENDENT_TERMS

RECURSIVE SUBROUTINE COMPUTE_TBA_TERM(Hamiltonian, kx, ky)
  IMPLICIT NONE
  REAL*8, INTENT(INOUT) :: kx, ky
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !Twice as big because of spin
  INTEGER*4 :: spin, lat, row, col
  !Only specifying upper triangle of matrix, since Hamiltonian is hermitian
  DO spin = 0, 1
    DO lat = 0, SUBLATTICES - 2
      !Kinetic hopping between nearest neighbours (different sublattices/ Ti layers)
      row = spin * TBA_DIM + lat * ORBITALS + 1
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 1

      Hamiltonian(row, col) = Hamiltonian(row, col) + epsilon_yz(kx, ky)
      Hamiltonian(row + 1, col + 1) = Hamiltonian(row + 1, col + 1) + epsilon_zx(kx, ky)
      Hamiltonian(row + 2, col + 2) = Hamiltonian(row + 2, col + 2) + epsilon_xy(kx, ky)
    END DO
  END DO

  !Nambu space: H(k) -> -H(-k)
  kx = -kx
  ky = -ky
  DO spin = 0, 1
    DO lat = 0, SUBLATTICES - 2
      !Kinetic hopping between nearest neighbours (different sublattices/ Ti layers)
      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 1

      Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(epsilon_yz(kx, ky))
      Hamiltonian(row + 1, col + 1) = Hamiltonian(row + 1, col + 1) - CONJG(epsilon_zx(kx, ky))
      Hamiltonian(row + 2, col + 2) = Hamiltonian(row + 2, col + 2) - CONJG(epsilon_xy(kx, ky))
    END DO
  END DO
  kx = -kx
  ky = -ky

END SUBROUTINE COMPUTE_TBA_TERM

RECURSIVE SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: nambu, lat, row, col
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO lat = 0, SUBLATTICES - 1
      row = nambu * DIM_POSITIVE_K + lat * ORBITALS + 1
      col = nambu * DIM_POSITIVE_K + lat * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) + imag * lambda_SOC / 2.

      row = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 1
      col = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - imag * lambda_SOC / 2.

      row = nambu * DIM_POSITIVE_K + lat * ORBITALS + 1
      col = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) - sign * lambda_SOC / 2. !Only real elements change sign here

      row = nambu * DIM_POSITIVE_K + lat * ORBITALS + 2
      col = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) + imag * lambda_SOC / 2.

      row = nambu * DIM_POSITIVE_K + lat * ORBITALS + 3
      col = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) + sign * lambda_SOC / 2. !Only real elements change sign here

      row = nambu * DIM_POSITIVE_K + lat * ORBITALS + 3
      col = nambu * DIM_POSITIVE_K + TBA_DIM + lat * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - imag * lambda_SOC / 2.
    END DO
  END DO

END SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS

RECURSIVE SUBROUTINE COMPUTE_TRIGONAL_TERMS(Hamiltonian)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: lat, spin, nambu, row, col
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO lat = 0, SUBLATTICES - 1
      DO spin = 0, 1
        row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
        col = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.

        row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
        col = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.

        row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
        col = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_TRIGONAL_TERMS

RECURSIVE SUBROUTINE COMPUTE_ELECTRIC_FIELD(Hamiltonian)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: i

  DO i = 1, 3
    !Ti1 atoms
    Hamiltonian(i, i) = Hamiltonian(i, i) + v / 2.
    Hamiltonian(i + TBA_DIM, i + TBA_DIM) = Hamiltonian(i + TBA_DIM, i + TBA_DIM) + v / 2.
    !Ti2 atoms
    Hamiltonian(i + ORBITALS, i + ORBITALS) = Hamiltonian(i + ORBITALS, i + ORBITALS) - v / 2.
    Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) = Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) - v / 2

    !Nambu space
    !Ti1 atoms
    Hamiltonian(DIM_POSITIVE_K + i, DIM_POSITIVE_K + i) = Hamiltonian(DIM_POSITIVE_K + i, DIM_POSITIVE_K + i) - v / 2.
    Hamiltonian(i + TBA_DIM + DIM_POSITIVE_K, i + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(i + TBA_DIM + DIM_POSITIVE_K, i + TBA_DIM + DIM_POSITIVE_K) - v / 2.
    !Ti2 atoms
    Hamiltonian(i + ORBITALS + DIM_POSITIVE_K, i + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(i + ORBITALS + DIM_POSITIVE_K, i + ORBITALS + DIM_POSITIVE_K) + v / 2.
    Hamiltonian(i + TBA_DIM + ORBITALS + DIM_POSITIVE_K, i + TBA_DIM + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(i + TBA_DIM + ORBITALS + DIM_POSITIVE_K, i + TBA_DIM + ORBITALS + DIM_POSITIVE_K) + v / 2

  END DO
END SUBROUTINE COMPUTE_ELECTRIC_FIELD

RECURSIVE SUBROUTINE COMPUTE_TI1_TI2(Hamiltonian, kx, ky)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  !Spin-up part
  Hamiltonian(1, 2 + ORBITALS) = Hamiltonian(1, 2 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
  Hamiltonian(1, 3 + ORBITALS) = Hamiltonian(1, 3 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
  Hamiltonian(2, 1 + ORBITALS) = Hamiltonian(2, 1 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
  Hamiltonian(2, 3 + ORBITALS) = Hamiltonian(2, 3 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))
  Hamiltonian(3, 1 + ORBITALS) = Hamiltonian(3, 1 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
  Hamiltonian(3, 2 + ORBITALS) = Hamiltonian(3, 2 + ORBITALS) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))

  !Spin-down part
  Hamiltonian(1 + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
  Hamiltonian(1 + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
  Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
  Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))
  Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
  Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))

  !Nambu space
  kx = -kx
  ky = -ky
  !Spin-up part
  Hamiltonian(1 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
  Hamiltonian(1 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
  Hamiltonian(2 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
  Hamiltonian(2 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
  Hamiltonian(3 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(3 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
  Hamiltonian(3 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(3 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))

  !Spin-down part
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
  Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
  Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
  Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
  Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.) * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
  kx = -kx
  ky = -ky

END SUBROUTINE COMPUTE_TI1_TI2

RECURSIVE SUBROUTINE COMPUTE_H_PI(Hamiltonian, kx, ky)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  REAL*8 :: k1, k2, k3

  k1 = -SQRT(3.) / 2.*kx + 3./2.*ky
  k2 = -SQRT(3.) / 2.*kx - 3./2.*ky
  k3 = SQRT(3.) * kx

  !Spin up, Ti1
  Hamiltonian(1, 2) = Hamiltonian(1, 2) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
  Hamiltonian(1, 3) = Hamiltonian(1, 3) + eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
  Hamiltonian(2, 3) = Hamiltonian(2, 3) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3))
  !Spin up, Ti2
  Hamiltonian(1 + ORBITALS, 2 + ORBITALS) = Hamiltonian(1 + ORBITALS, 2 + ORBITALS) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
  Hamiltonian(1 + ORBITALS, 3 + ORBITALS) = Hamiltonian(1 + ORBITALS, 3 + ORBITALS) + eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
  Hamiltonian(2 + ORBITALS, 3 + ORBITALS) = Hamiltonian(2 + ORBITALS, 3 + ORBITALS) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3))

  !Spin down, Ti1
  Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
  Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) + eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
  Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3))
  !Spin down, Ti2
  Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
  Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
  Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) - eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3))

  !Nambu space
  k1 = -k1
  k2 = -k2
  k3 = -k3
  !Spin up, Ti1
  Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
  Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
  Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3)))
  !Spin up, Ti2
  Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
  Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
  Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3)))

  !Spin down, Ti1
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
  Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3)))
  !Spin down, Ti2
  Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
  Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
  Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * V_pdp / SQRT(15.) * (2 * SIN(k1) + SIN(k2) + SIN(k3)))

END SUBROUTINE COMPUTE_H_PI

RECURSIVE SUBROUTINE COMPUTE_H_SIGMA(Hamiltonian, kx, ky)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  REAL*8 :: k1, k2, k3
  k1 = -SQRT(3.) / 2.*kx + 3./2.*ky
  k2 = -SQRT(3.) / 2.*kx - 3./2.*ky
  k3 = SQRT(3.) * kx

  !Spin up, Ti1
  Hamiltonian(1, 2) = Hamiltonian(1, 2) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(1, 3) = Hamiltonian(1, 3) - eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(2, 3) = Hamiltonian(2, 3) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3))
  !Spin up, Ti2
  Hamiltonian(1 + ORBITALS, 2 + ORBITALS) = Hamiltonian(1 + ORBITALS, 2 + ORBITALS) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(1 + ORBITALS, 3 + ORBITALS) = Hamiltonian(1 + ORBITALS, 3 + ORBITALS) - eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(2 + ORBITALS, 3 + ORBITALS) = Hamiltonian(2 + ORBITALS, 3 + ORBITALS) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3))
  !Spin down, Ti1
  Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) - eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3))
  !Spin down, Ti2
  Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) - eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2))
  Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3))

  !Nambu space
  k1 = -k1
  k2 = -k2
  k3 = -k3
  !Spin up, Ti1
  Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3)))
  !Spin up, Ti2
  Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3)))
  !Spin down, Ti1
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3)))
  !Spin down, Ti2
  Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k1) + SIN(k2)))
  Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.) * (SIN(k2) + SIN(k3)))

END SUBROUTINE COMPUTE_H_SIGMA

RECURSIVE SUBROUTINE COMPUTE_TETRAGONAL_STRAIN(Hamiltonian)
  !! Apply tetragonal strain effect. This energetically favours a single orbital.
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !! Hamiltonian that is to be modified

  INTEGER*4 :: nambu, spin, lat, row
  REAL*8 :: sign

  DO nambu = 0,1
    sign = (-1)**nambu
    DO spin = 0, 1
      DO lat = 0, SUBLATTICES - 1
        row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + orb_affected_tetragonal
        Hamiltonian(row, row) = Hamiltonian(row, row) + sign * zeta_tetragonal
      END DO
    END DO
  END DO
END SUBROUTINE COMPUTE_TETRAGONAL_STRAIN


SUBROUTINE COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  INTEGER*4 :: lat, spin, row, col

  DO spin = 0, 1
    DO lat = 0, SUBLATTICES - 2
      row = spin * TBA_DIM + lat * ORBITALS + 1
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_yz_zx(kx, ky)

      row = spin * TBA_DIM + lat * ORBITALS + 1
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_yz_xy(kx, ky)

      row = spin * TBA_DIM + lat * ORBITALS + 2
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_yz_zx(kx, ky)

      row = spin * TBA_DIM + lat * ORBITALS + 2
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_zx_xy(kx, ky)

      row = spin * TBA_DIM + lat * ORBITALS + 3
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_yz_xy(kx, ky)

      row = spin * TBA_DIM + lat * ORBITALS + 3
      col = spin * TBA_DIM + (lat + 1) * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_zx_xy(kx, ky)
    END DO
  END DO

  !Nambu space: H(k) -> -H(-k)
  kx = -kx
  ky = -ky
  DO spin = 0, 1
    DO lat = 0, SUBLATTICES - 2
      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_yz_zx(kx, ky))

      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_yz_xy(kx, ky))

      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_yz_zx(kx, ky))

      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_zx_xy(kx, ky))

      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_yz_xy(kx, ky))

      row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
      col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_zx_xy(kx, ky))
    END DO
  END DO
  kx = -kx
  ky = -ky
END SUBROUTINE COMPUTE_RASHBA_HOPPING

SUBROUTINE COMPUTE_LAYER_POTENTIAL(Hamiltonian)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: nambu, spin, lat, orb, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, 1
      DO lat = 0, SUBLATTICES - 1
        DO orb = 1, ORBITALS
          row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * V_layer(lat + 1) !lat + 1, because we cannot start from 0
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE COMPUTE_LAYER_POTENTIAL

SUBROUTINE COMPUTE_SUBBAND_POTENTIAL(Hamiltonian, n_band)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4, INTENT(IN) :: n_band
  INTEGER*4 :: nambu, row, i
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO i = 1, DIM_POSITIVE_K
      row = nambu * DIM_POSITIVE_K + i
      Hamiltonian(row, row) = Hamiltonian(row, row) + sign * Subband_energies(n_band) !lat + 1, because we cannot start from 0
    END DO
  END DO
END SUBROUTINE COMPUTE_SUBBAND_POTENTIAL

SUBROUTINE COMPUTE_FERMI_ENERGY(Hamiltonian)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  INTEGER*4 :: nambu, i, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO i = 1, DIM_POSITIVE_K
      row = nambu * DIM_POSITIVE_K + i
      Hamiltonian(row, row) = Hamiltonian(row, row) - sign * E_Fermi
    END DO
  END DO
END SUBROUTINE COMPUTE_FERMI_ENERGY

RECURSIVE SUBROUTINE COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC)
    !! Computes the superconducting coupling at given (kx,ky) point
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !! Hamiltonian of the system that is to be filled
  COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS) !! Superconducting energies
  REAL*8, INTENT(IN) :: kx !! Wavevector in X direction
  REAL*8, INTENT(IN) :: ky !! Wavevector in Y direction
  INTEGER*4 :: orb, lat, spin, row, col, gamma_spin_index, gamma_lat_index

  !Nearest neighbours pairing
  DO orb = 1, ORBITALS
    DO spin = 0, 1
      gamma_spin_index = MOD(spin + 1, 2) + 1
      ! -2, because we have to iterate up to one-before-last sublattice to be able to increment (lat + 1)
      DO lat = 0, SUBLATTICES - 2
        gamma_lat_index = 2 * lat + 2
        !Ti1 - Ti2 coupling
        row = spin * TBA_DIM + orb + lat * ORBITALS
        col = orb + (lat + 1) * ORBITALS + DIM_POSITIVE_K + TBA_DIM * MOD(spin + 1, 2)
        ! PRINT*, row, col
        Hamiltonian(row, col) = Hamiltonian(row, col) + &
        & Gamma_SC(orb, 1, gamma_spin_index, gamma_lat_index) * pairing_1(ky) +&
        & Gamma_SC(orb, 2, gamma_spin_index, gamma_lat_index) * pairing_2(kx, ky) +&
        & Gamma_SC(orb, 3, gamma_spin_index, gamma_lat_index) * pairing_3(kx, ky)

        gamma_lat_index = 2 * lat + 1

        !Ti2 - Ti1 coupling
        row = spin * TBA_DIM + orb + (lat + 1) * ORBITALS
        col = orb + DIM_POSITIVE_K + TBA_DIM * MOD(spin + 1, 2) + lat * ORBITALS
        ! PRINT*, row, col

        Hamiltonian(row, col) = Hamiltonian(row, col) + &
        & Gamma_SC(orb, 1, gamma_spin_index, gamma_lat_index) * CONJG(pairing_1(ky)) +&
        & Gamma_SC(orb, 2, gamma_spin_index, gamma_lat_index) * CONJG(pairing_2(kx, ky)) +&
        & Gamma_SC(orb, 3, gamma_spin_index, gamma_lat_index) * CONJG(pairing_3(kx, ky))
      END DO
    END DO
  END DO

  !Next nearest neighbours pairing
  DO orb = 1, ORBITALS
    DO lat = 0, SUBLATTICES - 1
      DO spin = 0, 1
        row = orb + lat * ORBITALS + spin * TBA_DIM
        col = orb + lat * ORBITALS + MOD(spin + 1, 2) * TBA_DIM + DIM_POSITIVE_K
        Hamiltonian(row, col) = Hamiltonian(row, col) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 1, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_1(kx)) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 2, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_2(kx, ky)) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 3, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_3(kx, ky)) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 4, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_4(kx)) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 5, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_5(kx, ky)) + &
        & Gamma_SC(orb, N_NEIGHBOURS + 6, MOD(spin + 1, 2) + 1, lat + 1) * CONJG(pairing_nnn_6(kx, ky))
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_SC

RECURSIVE SUBROUTINE COMPUTE_HUBBARD(Hamiltonian, Charge_dens)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K)
  INTEGER*4 :: orb, lat, orb_prime, spin, nambu, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, 1
      DO lat = 0, SUBLATTICES - 1
        DO orb = 1, ORBITALS
          row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * U_HUB * Charge_dens(MOD(spin + 1, 2) * TBA_DIM + lat * ORBITALS + orb) !Modulo should give opposite spin
          DO orb_prime = 1, ORBITALS
            IF (orb .NE. orb_prime) THEN
              Hamiltonian(row, row) = Hamiltonian(row, row) + sign * V_HUB * (Charge_dens(lat * ORBITALS + orb_prime) + Charge_dens(TBA_DIM + lat * ORBITALS + orb_prime)) !total charge dens in orbital
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO

  ! !Nambu space
  ! DO spin = 0, 1
  !     DO lat = 0, SUBLATTICES - 1
  !         DO orb = 1, ORBITALS
  !             Hamiltonian(DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb, DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb) = Hamiltonian(DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb, DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb) - &
  !             & U_HUB*Charge_dens(MOD(spin + 1, 2)*TBA_DIM + lat*ORBITALS + orb) !Modulo should give opposite spin
  !             DO orb_prime = 1, ORBITALS
  !                 IF (orb .NE. orb_prime) THEN
  !                     Hamiltonian(DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb, DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb) = Hamiltonian(DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb, DIM_POSITIVE_K + spin*TBA_DIM + lat*ORBITALS + orb) - &
  !                     & V_HUB*(Charge_dens(lat*ORBITALS + orb_prime) + Charge_dens(TBA_DIM + lat*ORBITALS + orb_prime)) !total charge dens in orbital
  !                 END IF
  !             END DO
  !         END DO
  !     END DO
  ! END DO

END SUBROUTINE COMPUTE_HUBBARD

RECURSIVE SUBROUTINE COMPUTE_ZEEMAN(B, Hamiltonian)
  REAL*8, INTENT(IN) :: B(3)
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, PARAMETER :: gFactor = 3.0d0
  REAL*8, PARAMETER :: muB = 0.5
  INTEGER*4 :: i, spin, nambu, row, col
  REAL*8 :: sign_nambu, sign_spin

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO spin = 0, 1
      sign_spin = (-1)**spin
      DO i = 1, TBA_DIM
        row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + i
        col = MIN(row + TBA_DIM, DIM)
        !B_z terms
        Hamiltonian(row, row) = Hamiltonian(row, row) + sign_nambu * sign_spin * 0.5d0 * muB * gFactor * B(3)
      END DO
    END DO
  END DO

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO i = 1, TBA_DIM
      row = nambu * DIM_POSITIVE_K + i
      col = nambu * DIM_POSITIVE_K + TBA_DIM + i
      !B_x and B_y terms
      Hamiltonian(row, col) = Hamiltonian(row, col) + sign_nambu * 0.5d0 * muB * gFactor * (B(1) - imag * B(2))
    END DO
  END DO
END SUBROUTINE COMPUTE_ZEEMAN

RECURSIVE SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING(B, Hamiltonian)
  !! Computes L \cdot B coupling, taking into account d orbitals
  REAL*8, INTENT(IN) :: B(3) !! Magnetic field
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM) !! Hamiltonian to be updated

  INTEGER*4 :: i_orb, j_orb, lat, spin, nambu, row, col
  REAL*8 :: sign_nambu
  REAL*8, PARAMETER :: muB = 0.5
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

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO spin = 0, 1
      DO lat = 0, SUBLATTICES - 1
        DO i_orb = 1, ORBITALS
          DO j_orb = i_orb, ORBITALS
            row = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + i_orb
            col = nambu * DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + j_orb
            Hamiltonian(row, col) = Hamiltonian(row, col) + sign_nambu * muB * &
            & (B(1) * L_x(i_orb, j_orb) + B(2) * L_y(i_orb, j_orb) + B(3) * L_z(i_orb, j_orb))
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING

!################ ADDITIONAL HAMILTONIANS FOR TESTING ########################
RECURSIVE SUBROUTINE DIAGONAL_TBA(Hamiltonian, kx, ky)
  IMPLICIT NONE
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM, DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  INTEGER*4 :: i

  DO i = 1, DIM_POSITIVE_K
    Hamiltonian(i, i) = Hamiltonian(i, i) + (DCOS(kx) + DCOS(ky))
  END DO

  !Nambu space
  kx = -kx
  ky = -ky
  DO i = DIM_POSITIVE_K + 1, DIM
    Hamiltonian(i, i) = Hamiltonian(i, i) - (DCOS(kx) + DCOS(ky))
  END DO
  kx = -kx
  ky = -ky

END SUBROUTINE

END MODULE mod_hamiltonians
