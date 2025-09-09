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

MODULE hamiltonians
USE utilities
USE parameters
USE reader
USE types
IMPLICIT NONE
CONTAINS

PURE RECURSIVE SUBROUTINE COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian, discretization, physical_params)
  !! Computes all terms that do not depend on k, including complex conjugate elements
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  TYPE(physical_params_t), INTENT(IN) :: physical_params
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
  CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian, physical_params % subband_params % delta_trigonal, discretization)
  CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, physical_params % subband_params % lambda_SOC, discretization)
  CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian, physical_params % subband_params % v, discretization)
  CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian, physical_params % subband_params % V_layer, discretization)
  CALL COMPUTE_TETRAGONAL_STRAIN(Hamiltonian, physical_params % subband_params % zeta_tetragonal, physical_params % subband_params % orb_affected_tetragonal, discretization)
  CALL COMPUTE_FERMI_ENERGY(Hamiltonian, physical_params % subband_params % E_Fermi, discretization)
  CALL COMPUTE_ZEEMAN(physical_params % external % B_field, physical_params % subband_params % g_factor, Hamiltonian, discretization)
  CALL COMPUTE_ORBITAL_MAGNETIC_COUPLING(physical_params % external % B_field, Hamiltonian, discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle
END SUBROUTINE COMPUTE_K_INDEPENDENT_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, discretization, physical_params)
  !! Computes all terms that depend on k, excluding complex conjugate elements
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  TYPE(physical_params_t), INTENT(IN) :: physical_params
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
  REAL*8, INTENT(INOUT) :: kx !! Wavevector in X direction
  REAL*8, INTENT(INOUT) :: ky !! Wavevector in Y direction
  CALL COMPUTE_TBA_TERM(Hamiltonian(:, :), kx, ky, physical_params % subband_params % t_D, physical_params % subband_params % t_I, discretization)
  CALL COMPUTE_TI1_TI2(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pdp, discretization)
  CALL COMPUTE_H_PI(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pdp, discretization)
  CALL COMPUTE_H_SIGMA(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pds, discretization)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
  CALL COMPUTE_RASHBA_HOPPING(Hamiltonian(:, :), kx, ky, physical_params % subband_params % t_Rashba, discretization) !This is adapted from KTaO_3, see: PRB, 103, 035115
END SUBROUTINE COMPUTE_K_DEPENDENT_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_TBA_TERM(Hamiltonian, kx, ky, t_D, t_I, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: t_D, t_I
  REAL*8, INTENT(IN) :: kx, ky
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, &
                                          & discretization % derived % DIM) !Twice as big because of spin
  INTEGER*4 :: spin, lat, row, col
  !Only specifying upper triangle of matrix, since Hamiltonian is hermitian
  DO spin = 0, SPINS - 1
    DO lat = 0, discretization % SUBLATTICES - 2
      !Kinetic hopping between nearest neighbours (different sublattices/ Ti layers)
      row = spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
      col = spin * discretization % derived % TBA_DIM + (lat + 1) * discretization % ORBITALS + 1

      Hamiltonian(row, col) = Hamiltonian(row, col) + epsilon_yz(kx, ky, t_D, t_I)
      Hamiltonian(row + 1, col + 1) = Hamiltonian(row + 1, col + 1) + epsilon_zx(kx, ky, t_D, t_I)
      Hamiltonian(row + 2, col + 2) = Hamiltonian(row + 2, col + 2) + epsilon_xy(kx, ky, t_D, t_I)
    END DO
  END DO

  !Nambu space: H(k) -> -H(-k)
  DO spin = 0, 1
    DO lat = 0, discretization % SUBLATTICES - 2
      !Kinetic hopping between nearest neighbours (different sublattices/ Ti layers)
      row = discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
      col = discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + (lat + 1) * discretization % ORBITALS + 1

      Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(epsilon_yz(-kx, -ky, t_D, t_I))
      Hamiltonian(row + 1, col + 1) = Hamiltonian(row + 1, col + 1) - CONJG(epsilon_zx(-kx, -ky, t_D, t_I))
      Hamiltonian(row + 2, col + 2) = Hamiltonian(row + 2, col + 2) - CONJG(epsilon_xy(-kx, -ky, t_D, t_I))
    END DO
  END DO

END SUBROUTINE COMPUTE_TBA_TERM

PURE RECURSIVE SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, lambda_SOC, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: lambda_SOC
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4 :: nambu, lat, row, col
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO lat = 0, discretization % SUBLATTICES - 1
      row = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 1
      col = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) + imag * lambda_SOC / 2.

      row = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - imag * lambda_SOC / 2.

      row = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 1
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) - sign * lambda_SOC / 2. !Only real elements change sign here

      row = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 2
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
      Hamiltonian(row, col) = Hamiltonian(row, col) + imag * lambda_SOC / 2.

      row = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 3
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
      Hamiltonian(row, col) = Hamiltonian(row, col) + sign * lambda_SOC / 2. !Only real elements change sign here

      row = nambu * discretization % derived % DIM_POSITIVE_K + lat * discretization % ORBITALS + 3
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
      Hamiltonian(row, col) = Hamiltonian(row, col) - imag * lambda_SOC / 2.
    END DO
  END DO

END SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_TRIGONAL_TERMS(Hamiltonian, delta_trigonal, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: delta_trigonal
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4 :: lat, spin, nambu, row, col
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO lat = 0, discretization % SUBLATTICES - 1
      DO spin = 0, 1
        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.

        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.

        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * delta_trigonal / 2.
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_TRIGONAL_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_ELECTRIC_FIELD(Hamiltonian, v, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: v
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4 :: nambu, orb, lat, spin
  INTEGER*4 :: row
  REAL*8 :: sign

  !! ATTENTION: This works only for two layers!!!!!
  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        sign = sign * (-1)**lat
        DO orb = 1, discretization % ORBITALS
          row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * v / 2.
        END DO
      END DO
    END DO
  END DO

  ! DO i = 1, 3
  !   !Ti1 atoms
  !   Hamiltonian(i, i) = Hamiltonian(i, i) + v / 2.
  !   Hamiltonian(i + TBA_DIM, i + TBA_DIM) = Hamiltonian(i + TBA_DIM, i + TBA_DIM) + v / 2.
  !   !Ti2 atoms
  !   Hamiltonian(i + ORBITALS, i + ORBITALS) = Hamiltonian(i + ORBITALS, i + ORBITALS) - v / 2.
  !   Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) = Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) - v / 2

  !   !Nambu space
  !   !Ti1 atoms
  !   Hamiltonian(DIM_POSITIVE_K + i, DIM_POSITIVE_K + i) = Hamiltonian(DIM_POSITIVE_K + i, DIM_POSITIVE_K + i) - v / 2.
  !   Hamiltonian(i + TBA_DIM + DIM_POSITIVE_K, i + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(i + TBA_DIM + DIM_POSITIVE_K, i + TBA_DIM + DIM_POSITIVE_K) - v / 2.
  !   !Ti2 atoms
  !   Hamiltonian(i + ORBITALS + DIM_POSITIVE_K, i + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(i + ORBITALS + DIM_POSITIVE_K, i + ORBITALS + DIM_POSITIVE_K) + v / 2.
  !   Hamiltonian(i + TBA_DIM + ORBITALS + DIM_POSITIVE_K, i + TBA_DIM + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(i + TBA_DIM + ORBITALS + DIM_POSITIVE_K, i + TBA_DIM + ORBITALS + DIM_POSITIVE_K) + v / 2

  ! END DO
END SUBROUTINE COMPUTE_ELECTRIC_FIELD

PURE RECURSIVE SUBROUTINE COMPUTE_TI1_TI2(Hamiltonian, kx, ky, eta_p, V_pdp, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: eta_p, V_pdp
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, INTENT(INOUT) :: kx, ky
  REAL*8 :: strength
  !TODO: Generalize this for many sublattices
  strength = eta_p * V_pdp * SQRT(2.)**(7./4.) / SQRT(15.)
  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM)
    !Spin-up part
    Hamiltonian(1, 2 + ORBITALS) = Hamiltonian(1, 2 + ORBITALS) + strength * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
    Hamiltonian(1, 3 + ORBITALS) = Hamiltonian(1, 3 + ORBITALS) + strength * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
    Hamiltonian(2, 1 + ORBITALS) = Hamiltonian(2, 1 + ORBITALS) + strength * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
    Hamiltonian(2, 3 + ORBITALS) = Hamiltonian(2, 3 + ORBITALS) + strength * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))
    Hamiltonian(3, 1 + ORBITALS) = Hamiltonian(3, 1 + ORBITALS) + strength * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
    Hamiltonian(3, 2 + ORBITALS) = Hamiltonian(3, 2 + ORBITALS) + strength * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))

    !Spin-down part
    Hamiltonian(1 + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + strength * (-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
    Hamiltonian(1 + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + strength * (1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
    Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + strength * (2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx))
    Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + strength * (1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))
    Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + strength * (-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky)))
    Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + strength * (-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky)))

    !Nambu space
    kx = -kx
    ky = -ky
    !Spin-up part
    Hamiltonian(1 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
    Hamiltonian(1 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
    Hamiltonian(2 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
    Hamiltonian(2 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
    Hamiltonian(3 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(3 + DIM_POSITIVE_K, 1 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
    Hamiltonian(3 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(3 + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - strength * CONJG((-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))

    !Spin-down part
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((-2.*imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((1 - EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
    Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((2 * imag * EXP(imag * 3./2.*ky) * SIN(SQRT(3.) / 2.*kx)))
    Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((1 - EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
    Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((-1 + EXP(imag / 2.*(SQRT(3.) * kx + 3.*ky))))
    Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(3 + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - strength * CONJG((-1 + EXP(-imag / 2.*(SQRT(3.) * kx - 3.*ky))))
    kx = -kx
    ky = -ky
  END ASSOCIATE
END SUBROUTINE COMPUTE_TI1_TI2

PURE RECURSIVE SUBROUTINE COMPUTE_H_PI(Hamiltonian, kx, ky, eta_p, V_pdp, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: eta_p, V_pdp
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, INTENT(IN) :: kx, ky
  REAL*8 :: k1, k2, k3
  COMPLEX*16 :: strength

  !TODO: Generalize this for many sublattices
  strength = eta_p * 2 * imag * V_pdp / SQRT(15.)

  k1 = -SQRT(3.) / 2.*kx + 3./2.*ky
  k2 = -SQRT(3.) / 2.*kx - 3./2.*ky
  k3 = SQRT(3.) * kx
  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
            & ORBITALS => discretization % ORBITALS, &
            & TBA_DIM => discretization % derived % TBA_DIM)
    !Spin up, Ti1
    Hamiltonian(1, 2) = Hamiltonian(1, 2) - strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
    Hamiltonian(1, 3) = Hamiltonian(1, 3) + strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
    Hamiltonian(2, 3) = Hamiltonian(2, 3) - strength * (2 * SIN(k1) + SIN(k2) + SIN(k3))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS, 2 + ORBITALS) = Hamiltonian(1 + ORBITALS, 2 + ORBITALS) - strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
    Hamiltonian(1 + ORBITALS, 3 + ORBITALS) = Hamiltonian(1 + ORBITALS, 3 + ORBITALS) + strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
    Hamiltonian(2 + ORBITALS, 3 + ORBITALS) = Hamiltonian(2 + ORBITALS, 3 + ORBITALS) - strength * (2 * SIN(k1) + SIN(k2) + SIN(k3))

    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) - strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
    Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) + strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
    Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) - strength * (2 * SIN(k1) + SIN(k2) + SIN(k3))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) - strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3))
    Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) + strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3))
    Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) - strength * (2 * SIN(k1) + SIN(k2) + SIN(k3))

    !Nambu space
    k1 = -k1
    k2 = -k2
    k3 = -k3
    !Spin up, Ti1
    Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
    Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
    Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) + CONJG(strength * (2 * SIN(k1) + SIN(k2) + SIN(k3)))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
    Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
    Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) + CONJG(strength * (2 * SIN(k1) + SIN(k2) + SIN(k3)))

    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
    Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (2 * SIN(k1) + SIN(k2) + SIN(k3)))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2) + 2 * SIN(k3)))
    Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + 2 * SIN(k2) + SIN(k3)))
    Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (2 * SIN(k1) + SIN(k2) + SIN(k3)))
  END ASSOCIATE
END SUBROUTINE COMPUTE_H_PI

PURE RECURSIVE SUBROUTINE COMPUTE_H_SIGMA(Hamiltonian, kx, ky, eta_p, V_pds, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: eta_p, V_pds
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, INTENT(IN) :: kx, ky
  REAL*8 :: k1, k2, k3
  COMPLEX*16 :: strength

  !TODO: Generalize ths for many sublattices

  strength = eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.)

  k1 = -SQRT(3.) / 2.*kx + 3./2.*ky
  k2 = -SQRT(3.) / 2.*kx - 3./2.*ky
  k3 = SQRT(3.) * kx
  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM)
    !Spin up, Ti1
    Hamiltonian(1, 2) = Hamiltonian(1, 2) + strength * (SIN(k1) + SIN(k2))
    Hamiltonian(1, 3) = Hamiltonian(1, 3) - strength * (SIN(k1) + SIN(k2))
    Hamiltonian(2, 3) = Hamiltonian(2, 3) + strength * (SIN(k2) + SIN(k3))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS, 2 + ORBITALS) = Hamiltonian(1 + ORBITALS, 2 + ORBITALS) + strength * (SIN(k1) + SIN(k2))
    Hamiltonian(1 + ORBITALS, 3 + ORBITALS) = Hamiltonian(1 + ORBITALS, 3 + ORBITALS) - strength * (SIN(k1) + SIN(k2))
    Hamiltonian(2 + ORBITALS, 3 + ORBITALS) = Hamiltonian(2 + ORBITALS, 3 + ORBITALS) + strength * (SIN(k2) + SIN(k3))
    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + TBA_DIM) + strength * (SIN(k1) + SIN(k2))
    Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + TBA_DIM) - strength * (SIN(k1) + SIN(k2))
    Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + TBA_DIM) + strength * (SIN(k2) + SIN(k3))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 2 + ORBITALS + TBA_DIM) + strength * (SIN(k1) + SIN(k2))
    Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) - strength * (SIN(k1) + SIN(k2))
    Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM, 3 + ORBITALS + TBA_DIM) + strength * (SIN(k2) + SIN(k3))

    !Nambu space
    k1 = -k1
    k2 = -k2
    k3 = -k3
    !Spin up, Ti1
    Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 2 + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(1 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) = Hamiltonian(2 + DIM_POSITIVE_K, 3 + DIM_POSITIVE_K) - CONJG(strength * (SIN(k2) + SIN(k3)))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 2 + ORBITALS + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + DIM_POSITIVE_K, 3 + ORBITALS + DIM_POSITIVE_K) - CONJG(strength * (SIN(k2) + SIN(k3)))
    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 2 + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + TBA_DIM + DIM_POSITIVE_K, 3 + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k2) + SIN(k3)))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(1 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) + CONJG(strength * (SIN(k1) + SIN(k2)))
    Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) = Hamiltonian(2 + ORBITALS + TBA_DIM + DIM_POSITIVE_K, 3 + ORBITALS + TBA_DIM + DIM_POSITIVE_K) - CONJG(strength * (SIN(k2) + SIN(k3)))
  END ASSOCIATE
END SUBROUTINE COMPUTE_H_SIGMA

PURE RECURSIVE SUBROUTINE COMPUTE_TETRAGONAL_STRAIN(Hamiltonian, zeta_tetragonal, orb_affected_tetragonal, discretization)
  !! Apply tetragonal strain effect. This energetically favours a single orbital.
  IMPLICIT NONE

  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: zeta_tetragonal
  INTEGER*4, INTENT(IN) :: orb_affected_tetragonal
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian that is to be modified

  INTEGER*4 :: nambu, spin, lat, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb_affected_tetragonal
        Hamiltonian(row, row) = Hamiltonian(row, row) + sign * zeta_tetragonal
      END DO
    END DO
  END DO
END SUBROUTINE COMPUTE_TETRAGONAL_STRAIN

PURE SUBROUTINE COMPUTE_RASHBA_HOPPING(Hamiltonian, kx, ky, t_Rashba, discretization)
  IMPLICIT NONE

  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: t_Rashba
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, INTENT(IN) :: kx, ky
  INTEGER*4 :: lat, spin, row, col

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
            & TBA_DIM => discretization % derived % TBA_DIM, &
            & ORBITALS => discretization % ORBITALS, &
            & SUBLATTICES => discretization % SUBLATTICES)
    DO spin = 0, SPINS - 1
      DO lat = 0, SUBLATTICES - 2
        row = spin * TBA_DIM + lat * ORBITALS + 1
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_yz_zx(kx, ky, t_Rashba)

        row = spin * TBA_DIM + lat * ORBITALS + 1
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_yz_xy(kx, ky, t_Rashba)

        row = spin * TBA_DIM + lat * ORBITALS + 2
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 1
        Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_yz_zx(kx, ky, t_Rashba)

        row = spin * TBA_DIM + lat * ORBITALS + 2
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + rashba_zx_xy(kx, ky, t_Rashba)

        row = spin * TBA_DIM + lat * ORBITALS + 3
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 1
        Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_yz_xy(kx, ky, t_Rashba)

        row = spin * TBA_DIM + lat * ORBITALS + 3
        col = spin * TBA_DIM + (lat + 1) * ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) - rashba_zx_xy(kx, ky, t_Rashba)
      END DO
    END DO

    !Nambu space: H(k) -> -H(-k)
    DO spin = 0, 1
      DO lat = 0, SUBLATTICES - 2
        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))

        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 1
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))

        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 1
        Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_yz_zx(-kx, -ky, t_Rashba))

        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 2
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) - CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))

        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 1
        Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_yz_xy(-kx, -ky, t_Rashba))

        row = DIM_POSITIVE_K + spin * TBA_DIM + lat * ORBITALS + 3
        col = DIM_POSITIVE_K + spin * TBA_DIM + (lat + 1) * ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) + CONJG(rashba_zx_xy(-kx, -ky, t_Rashba))
      END DO
    END DO
  END ASSOCIATE
END SUBROUTINE COMPUTE_RASHBA_HOPPING

PURE SUBROUTINE COMPUTE_LAYER_POTENTIAL(Hamiltonian, V_layer, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: V_layer(discretization % SUBLATTICES)
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4 :: nambu, spin, lat, orb, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        DO orb = 1, discretization % ORBITALS
          row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * V_layer(lat + 1) !lat + 1, because we cannot start from 0
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE COMPUTE_LAYER_POTENTIAL

PURE SUBROUTINE COMPUTE_SUBBAND_POTENTIAL(Hamiltonian, n_band, Subband_energies, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: Subband_energies(discretization % SUBBANDS)
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4, INTENT(IN) :: n_band
  INTEGER*4 :: nambu, row, i
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO i = 1, discretization % derived % DIM_POSITIVE_K
      row = nambu * discretization % derived % DIM_POSITIVE_K + i
      Hamiltonian(row, row) = Hamiltonian(row, row) + sign * Subband_energies(n_band) !lat + 1, because we cannot start from 0
    END DO
  END DO
END SUBROUTINE COMPUTE_SUBBAND_POTENTIAL

PURE SUBROUTINE COMPUTE_FERMI_ENERGY(Hamiltonian, E_Fermi, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: E_Fermi
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER*4 :: nambu, i, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO i = 1, discretization % derived % DIM_POSITIVE_K
      row = nambu * discretization % derived % DIM_POSITIVE_K + i
      Hamiltonian(row, row) = Hamiltonian(row, row) - sign * E_Fermi
    END DO
  END DO
END SUBROUTINE COMPUTE_FERMI_ENERGY

PURE RECURSIVE SUBROUTINE COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, discretization)
    !! Computes the superconducting coupling at given (kx,ky) point
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
  COMPLEX*16, INTENT(IN) :: Gamma_SC(discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, discretization % derived % LAYER_COUPLINGS) !! Superconducting energies
  REAL*8, INTENT(IN) :: kx !! Wavevector in X direction
  REAL*8, INTENT(IN) :: ky !! Wavevector in Y direction
  INTEGER*4 :: orb, lat, spin1, spin2, row, col, gamma_lat_index
  !TODO: Change it so that it can use triplet couplings
  !Nearest neighbours pairing
  DO orb = 1, discretization % ORBITALS
    DO spin1 = 0, SPINS - 1
      DO spin2 = 0, SPINS - 1
        ! -2, because we have to iterate up to one-before-last sublattice to be able to increment (lat + 1)
        DO lat = 0, discretization % SUBLATTICES - 2
          gamma_lat_index = 2 * lat + 2
          !Ti1 - Ti2 coupling
          row = spin1 * discretization % derived % TBA_DIM + orb + lat * discretization % ORBITALS
          col = orb + (lat + 1) * discretization % ORBITALS + discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM * spin2
          ! PRINT*, row, col
          Hamiltonian(row, col) = Hamiltonian(row, col) + &
          & Gamma_SC(orb, 1, spin1 + 1, spin2 + 1, gamma_lat_index) * CONJG(pairing_1(ky)) +&
          & Gamma_SC(orb, 2, spin1 + 1, spin2 + 1, gamma_lat_index) * CONJG(pairing_2(kx, ky)) +&
          & Gamma_SC(orb, 3, spin1 + 1, spin2 + 1, gamma_lat_index) * CONJG(pairing_3(kx, ky))

          gamma_lat_index = 2 * lat + 1

          !Ti2 - Ti1 coupling
          row = spin1 * discretization % derived % TBA_DIM + orb + (lat + 1) * discretization % ORBITALS
          col = orb + discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM * spin2 + lat * discretization % ORBITALS
          ! PRINT*, row, col

          Hamiltonian(row, col) = Hamiltonian(row, col) + &
          & Gamma_SC(orb, 1, spin1 + 1, spin2 + 1, gamma_lat_index) * pairing_1(ky) +&
          & Gamma_SC(orb, 2, spin1 + 1, spin2 + 1, gamma_lat_index) * pairing_2(kx, ky) +&
          & Gamma_SC(orb, 3, spin1 + 1, spin2 + 1, gamma_lat_index) * pairing_3(kx, ky)
        END DO
      END DO
    END DO
  END DO
  !Next nearest neighbours pairing
  DO orb = 1, discretization % ORBITALS
    DO lat = 0, discretization % SUBLATTICES - 1
      DO spin1 = 0, SPINS - 1
        DO spin2 = 0, SPINS - 1
          row = orb + lat * discretization % ORBITALS + spin1 * discretization % derived % TBA_DIM
          col = orb + lat * discretization % ORBITALS + spin2 * discretization % derived % TBA_DIM + discretization % derived % DIM_POSITIVE_K
          Hamiltonian(row, col) = Hamiltonian(row, col) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 1, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_1(kx) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 2, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_2(kx, ky) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 3, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_3(kx, ky) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 4, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_4(kx) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 5, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_5(kx, ky) + &
          & Gamma_SC(orb, N_NEIGHBOURS + 6, spin1 + 1, spin2 + 1, lat + 1) * pairing_nnn_6(kx, ky)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_SC

PURE RECURSIVE SUBROUTINE COMPUTE_HUBBARD(Hamiltonian, Charge_dens, U_HUB, V_HUB, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: U_HUB, V_HUB
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K)
  INTEGER*4 :: orb, lat, orb_prime, spin, nambu, row
  REAL*8 :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        DO orb = 1, discretization % ORBITALS
          row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * U_HUB * Charge_dens(MOD(spin + 1, 2) * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb) !Modulo should give opposite spin
          DO orb_prime = 1, discretization % ORBITALS
            IF (orb .NE. orb_prime) THEN
              Hamiltonian(row, row) = Hamiltonian(row, row) + sign * V_HUB * (Charge_dens(lat * discretization % ORBITALS + orb_prime) + &
                & Charge_dens(discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb_prime)) !total charge dens in orbital
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_HUBBARD

PURE RECURSIVE SUBROUTINE COMPUTE_ZEEMAN(B, g_factor, Hamiltonian, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: B(3)
  REAL*8, INTENT(IN) :: g_factor
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL*8, PARAMETER :: muB = 0.5
  INTEGER*4 :: i, spin, nambu, row, col
  REAL*8 :: sign_nambu, sign_spin

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO spin = 0, SPINS - 1
      sign_spin = (-1)**spin
      DO i = 1, discretization % derived % TBA_DIM
        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + i
        col = MIN(row + discretization % derived % TBA_DIM, discretization % derived % DIM)
        !B_z terms
        Hamiltonian(row, row) = Hamiltonian(row, row) + sign_nambu * sign_spin * 0.5d0 * muB * g_factor * B(3)
      END DO
    END DO
  END DO

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO i = 1, discretization % derived % TBA_DIM
      row = nambu * discretization % derived % DIM_POSITIVE_K + i
      col = nambu * discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + i
      !B_x and B_y terms
      Hamiltonian(row, col) = Hamiltonian(row, col) + sign_nambu * 0.5d0 * muB * g_factor * (B(1) - imag * B(2))
    END DO
  END DO
END SUBROUTINE COMPUTE_ZEEMAN

PURE RECURSIVE SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING(B, Hamiltonian, discretization)
  !! Computes L \cdot B coupling, taking into account d orbitals
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: B(3) !! Magnetic field
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian to be updated

  INTEGER*4 :: i_orb, j_orb, lat, spin, nambu, row, col
  REAL*8 :: sign_nambu
  REAL*8, PARAMETER :: muB = 0.5
  COMPLEX*16, PARAMETER :: c_zero = (0.0d0, 0.0d0) ! So that compiler does not complain about type mismatch between 0 and imag
  COMPLEX*16, PARAMETER :: L_x(3, 3) = TRANSPOSE(RESHAPE([c_zero, c_zero, imag, &
                                                          c_zero, c_zero, imag, &
                                                          -imag, -imag, c_zero], &
                                                         [3, 3])) / SQRT(2.0d0)
  COMPLEX*16, PARAMETER :: L_y(3, 3) = TRANSPOSE(RESHAPE([c_zero, -2 * imag, -imag, &
                                                          2 * imag, c_zero, imag, &
                                                          imag, -imag, c_zero], &
                                                         [3, 3])) / SQRT(6.0d0)
  COMPLEX*16, PARAMETER :: L_z(3, 3) = TRANSPOSE(RESHAPE([c_zero, -imag, imag, &
                                                          imag, c_zero, -imag, &
                                                          -imag, imag, c_zero], &
                                                         [3, 3])) / SQRT(3.0d0)

  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        DO i_orb = 1, discretization % ORBITALS
          DO j_orb = i_orb, discretization % ORBITALS
            row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + i_orb
            col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + j_orb
            Hamiltonian(row, col) = Hamiltonian(row, col) + sign_nambu * muB * &
            & (B(1) * L_x(i_orb, j_orb) + B(2) * L_y(i_orb, j_orb) + B(3) * L_z(i_orb, j_orb))
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING

END MODULE hamiltonians
