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
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
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
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
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
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
  REAL(REAL64), INTENT(INOUT) :: kx !! Wavevector in X direction
  REAL(REAL64), INTENT(INOUT) :: ky !! Wavevector in Y direction
  CALL COMPUTE_TBA_TERM(Hamiltonian(:, :), kx, ky, physical_params % subband_params % t_D, physical_params % subband_params % t_I, discretization)
  CALL COMPUTE_TI1_TI2(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pdp, discretization)
  CALL COMPUTE_H_PI(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pdp, discretization)
  CALL COMPUTE_H_SIGMA(Hamiltonian(:, :), kx, ky, physical_params % subband_params % eta_p, physical_params % subband_params % V_pds, discretization)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
  CALL COMPUTE_RASHBA_HOPPING(Hamiltonian(:, :), kx, ky, physical_params % subband_params % t_Rashba, discretization) !This is adapted from KTaO_3, see: PRB, 103, 035115
END SUBROUTINE COMPUTE_K_DEPENDENT_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_TBA_TERM(Hamiltonian, kx, ky, t_D, t_I, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: t_D, t_I
  REAL(REAL64), INTENT(IN) :: kx, ky
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, &
                                          & discretization % derived % DIM) !Twice as big because of spin
  ! INTEGER(INT32) :: spin, lat, row, col

  COMPLEX(REAL64) :: Hoppings(discretization % ORBITALS * 2)

  INTEGER(INT32) :: nambu, spin, lat, current_lat_block, next_lat_block, spin_block, nambu_block, hopping_start_idx

  CALL COMPUTE_NEAREST_EVEN_HOPPING(Hoppings, kx, ky, t_D, t_I, discretization % ORBITALS * 2)

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM, &
             SUBLATTICES => discretization % SUBLATTICES)
    DO nambu = 0, 1
      nambu_block = nambu * DIM_POSITIVE_K
      hopping_start_idx = nambu * ORBITALS
      DO spin = 0, SPINS - 1
        spin_block = spin * TBA_DIM
        DO lat = 0, SUBLATTICES - 2
          current_lat_block = lat * ORBITALS + spin_block + nambu_block
          next_lat_block = current_lat_block + ORBITALS
          !Loop over orbitals unrolled manually
          Hamiltonian(current_lat_block + 1, next_lat_block + 1) = Hamiltonian(current_lat_block + 1, next_lat_block + 1) + Hoppings(hopping_start_idx + 1)
          Hamiltonian(current_lat_block + 2, next_lat_block + 2) = Hamiltonian(current_lat_block + 2, next_lat_block + 2) + Hoppings(hopping_start_idx + 2)
          Hamiltonian(current_lat_block + 3, next_lat_block + 3) = Hamiltonian(current_lat_block + 3, next_lat_block + 3) + Hoppings(hopping_start_idx + 3)
        END DO
      END DO
    END DO
  END ASSOCIATE

END SUBROUTINE COMPUTE_TBA_TERM

PURE RECURSIVE SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, lambda_SOC, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: lambda_SOC
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: nambu, lat, row, col, i_orb, j_orb, i_spin, j_spin
  REAL(REAL64) :: sign

  !H(k)
  DO i_spin = 1, SPINS
    DO j_spin = 1, SPINS
      DO lat = 0, discretization % SUBLATTICES - 1
        DO i_orb = 1, discretization % ORBITALS
          DO j_orb = 1, discretization % ORBITALS
            row = (i_spin - 1) * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + i_orb
            col = (j_spin - 1) * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + j_orb
            Hamiltonian(row, col) = Hamiltonian(row, col) - 0.5 * lambda_SOC * &
            & (L_x(i_orb, j_orb) * Sigma_x(i_spin, j_spin) + L_y(i_orb, j_orb) * Sigma_y(i_spin, j_spin) + L_z(i_orb, j_orb) * Sigma_z(i_spin, j_spin))
          END DO
        END DO
      END DO
    END DO
  END DO
  !-H*(-k)
  DO i_spin = 1, SPINS
    DO j_spin = 1, SPINS
      DO lat = 0, discretization % SUBLATTICES - 1
        DO i_orb = 1, discretization % ORBITALS
          DO j_orb = 1, discretization % ORBITALS
            row = discretization % derived % DIM_POSITIVE_K + (i_spin - 1) * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + i_orb
            col = discretization % derived % DIM_POSITIVE_K + (j_spin - 1) * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + j_orb
            Hamiltonian(row, col) = Hamiltonian(row, col) + 0.5 * lambda_SOC * &
            & CONJG(L_x(i_orb, j_orb) * Sigma_x(i_spin, j_spin) + L_y(i_orb, j_orb) * Sigma_y(i_spin, j_spin) + L_z(i_orb, j_orb) * Sigma_z(i_spin, j_spin))
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_TRIGONAL_TERMS(Hamiltonian, delta_trigonal, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: delta_trigonal
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: lat, spin, nambu, row, col
  REAL(REAL64) :: sign

  DO nambu = 0, 1
    sign = (-1)**nambu
    DO lat = 0, discretization % SUBLATTICES - 1
      DO spin = 0, 1
        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * 0.5d0 * delta_trigonal

        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 1
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * 0.5d0 * delta_trigonal

        row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 2
        col = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + 3
        Hamiltonian(row, col) = Hamiltonian(row, col) + sign * 0.5d0 * delta_trigonal
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_TRIGONAL_TERMS

PURE RECURSIVE SUBROUTINE COMPUTE_ELECTRIC_FIELD(Hamiltonian, v, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: v
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: nambu, orb, lat, spin
  INTEGER(INT32) :: row
  REAL(REAL64) :: sign, sign_nambu

  !! ATTENTION: This works only for two layers!!!!!
  DO nambu = 0, 1
    sign_nambu = (-1)**nambu
    DO spin = 0, SPINS - 1
      DO lat = 0, discretization % SUBLATTICES - 1
        sign = sign_nambu * (-1)**lat
        DO orb = 1, discretization % ORBITALS
          row = nambu * discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + orb
          Hamiltonian(row, row) = Hamiltonian(row, row) + sign * v / 2.
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE COMPUTE_ELECTRIC_FIELD

PURE RECURSIVE SUBROUTINE COMPUTE_TI1_TI2(Hamiltonian, kx, ky, eta_p, V_pdp, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: eta_p, V_pdp
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL(REAL64), INTENT(INOUT) :: kx, ky
  REAL(REAL64) :: strength
  COMPLEX(REAL64) :: Hoppings(discretization % ORBITALS * 2)

  INTEGER(INT32) :: nambu, spin, lat, current_lat_block, next_lat_block, spin_block, nambu_block, hopping_start_idx

  strength = eta_p * V_pdp * SQRT(2.)**(0.5 * 7.) / SQRT(15.)

  CALL COMPUTE_NEAREST_ODD_HOPPING(Hoppings, kx, ky, strength, discretization % ORBITALS * 2)

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM, &
             SUBLATTICES => discretization % SUBLATTICES)
    DO nambu = 0, 1
      nambu_block = nambu * DIM_POSITIVE_K
      hopping_start_idx = nambu * ORBITALS
      DO spin = 0, SPINS - 1
        spin_block = spin * TBA_DIM
        DO lat = 0, SUBLATTICES - 2
          current_lat_block = lat * ORBITALS + spin_block + nambu_block
          next_lat_block = current_lat_block + ORBITALS

          !Manually unrolling the double loop over orbitals - the number of orbitals is small and hardcoded anyways.
          Hamiltonian(1 + current_lat_block, 2 + next_lat_block) = Hamiltonian(1 + current_lat_block, 2 + next_lat_block) + Hoppings(hopping_start_idx + 1)
          Hamiltonian(1 + current_lat_block, 3 + next_lat_block) = Hamiltonian(1 + current_lat_block, 3 + next_lat_block) + Hoppings(hopping_start_idx + 2)
          Hamiltonian(2 + current_lat_block, 3 + next_lat_block) = Hamiltonian(2 + current_lat_block, 3 + next_lat_block) + Hoppings(hopping_start_idx + 3)
          Hamiltonian(2 + current_lat_block, 1 + next_lat_block) = Hamiltonian(2 + current_lat_block, 1 + next_lat_block) - Hoppings(hopping_start_idx + 1)
          Hamiltonian(3 + current_lat_block, 1 + next_lat_block) = Hamiltonian(3 + current_lat_block, 1 + next_lat_block) - Hoppings(hopping_start_idx + 2)
          Hamiltonian(3 + current_lat_block, 2 + next_lat_block) = Hamiltonian(3 + current_lat_block, 2 + next_lat_block) - Hoppings(hopping_start_idx + 3)

        END DO
      END DO
    END DO
  END ASSOCIATE
END SUBROUTINE COMPUTE_TI1_TI2

PURE RECURSIVE SUBROUTINE COMPUTE_H_PI(Hamiltonian, kx, ky, eta_p, V_pdp, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: eta_p, V_pdp
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: kx, ky
  COMPLEX(REAL64) :: strength
  COMPLEX(REAL64) :: Hoppings(discretization % ORBITALS * 2)

  INTEGER(INT32) :: nambu, spin, lat, lat_block, spin_block, nambu_block, hopping_start_idx

  strength = eta_p * 2 * imag * V_pdp / SQRT(15.)

  CALL COMPUTE_NEXT_PI_ODD_HOPPING(Hoppings, kx, ky, strength, discretization % ORBITALS * 2)

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
            & ORBITALS => discretization % ORBITALS, &
            & TBA_DIM => discretization % derived % TBA_DIM, &
            & SUBLATTICES => discretization % SUBLATTICES)

    DO nambu = 0, 1
      nambu_block = nambu * DIM_POSITIVE_K
      hopping_start_idx = nambu * ORBITALS
      DO spin = 0, SPINS - 1
        spin_block = spin * TBA_DIM
        DO lat = 0, SUBLATTICES - 1
          lat_block = lat * ORBITALS + spin_block + nambu_block

          !Manually unrolling the loop over orbitals - the number of orbitals is small and hardcoded anyways.
          Hamiltonian(1 + lat_block, 2 + lat_block) = Hamiltonian(1 + lat_block, 2 + lat_block) + Hoppings(hopping_start_idx + 1)
          Hamiltonian(1 + lat_block, 3 + lat_block) = Hamiltonian(1 + lat_block, 3 + lat_block) + Hoppings(hopping_start_idx + 2)
          Hamiltonian(2 + lat_block, 3 + lat_block) = Hamiltonian(2 + lat_block, 3 + lat_block) + Hoppings(hopping_start_idx + 3)
        END DO
      END DO
    END DO
  END ASSOCIATE
END SUBROUTINE COMPUTE_H_PI

PURE RECURSIVE SUBROUTINE COMPUTE_H_SIGMA(Hamiltonian, kx, ky, eta_p, V_pds, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: eta_p, V_pds
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: kx, ky
  COMPLEX(REAL64) :: strength
  COMPLEX(REAL64) :: Hoppings(discretization % ORBITALS * 2)

  INTEGER(INT32) :: nambu, spin, lat, lat_block, spin_block, nambu_block, hopping_start_idx

  strength = eta_p * 2 * imag * SQRT(3.) * V_pds / SQRT(15.)

  CALL COMPUTE_NEXT_SIGMA_ODD_HOPPING(Hoppings, kx, ky, strength, discretization % ORBITALS * 2)

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM, &
             SUBLATTICES => discretization % SUBLATTICES)

    DO nambu = 0, 1
      nambu_block = nambu * DIM_POSITIVE_K
      hopping_start_idx = nambu * ORBITALS
      DO spin = 0, SPINS - 1
        spin_block = spin * TBA_DIM
        DO lat = 0, SUBLATTICES - 1
          lat_block = lat * ORBITALS + spin_block + nambu_block

          !Manually unrolling the loop over orbitals - the number of orbitals is small and hardcoded anyways.
          Hamiltonian(1 + lat_block, 2 + lat_block) = Hamiltonian(1 + lat_block, 2 + lat_block) + Hoppings(hopping_start_idx + 1)
          Hamiltonian(1 + lat_block, 3 + lat_block) = Hamiltonian(1 + lat_block, 3 + lat_block) + Hoppings(hopping_start_idx + 2)
          Hamiltonian(2 + lat_block, 3 + lat_block) = Hamiltonian(2 + lat_block, 3 + lat_block) + Hoppings(hopping_start_idx + 3)
        END DO
      END DO
    END DO

  END ASSOCIATE
END SUBROUTINE COMPUTE_H_SIGMA

PURE RECURSIVE SUBROUTINE COMPUTE_TETRAGONAL_STRAIN(Hamiltonian, zeta_tetragonal, orb_affected_tetragonal, discretization)
  !! Apply tetragonal strain effect. This energetically favours a single orbital.
  IMPLICIT NONE

  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: zeta_tetragonal
  INTEGER(INT32), INTENT(IN) :: orb_affected_tetragonal
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian that is to be modified

  INTEGER(INT32) :: nambu, spin, lat, row
  REAL(REAL64) :: sign

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
  REAL(REAL64), INTENT(IN) :: t_Rashba
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL(REAL64), INTENT(INOUT) :: kx, ky

  COMPLEX(REAL64) :: Hoppings(discretization % ORBITALS * 2)

  INTEGER(INT32) :: nambu, spin, lat, current_lat_block, next_lat_block, spin_block, nambu_block, hopping_start_idx

  CALL COMPUTE_NEAREST_ODD_HOPPING(Hoppings, kx, ky, t_Rashba, discretization % ORBITALS * 2)

  ASSOCIATE (DIM_POSITIVE_K => discretization % derived % DIM_POSITIVE_K, &
             ORBITALS => discretization % ORBITALS, &
             TBA_DIM => discretization % derived % TBA_DIM, &
             SUBLATTICES => discretization % SUBLATTICES)
    DO nambu = 0, 1
      nambu_block = nambu * DIM_POSITIVE_K
      hopping_start_idx = nambu * ORBITALS
      DO spin = 0, SPINS - 1
        spin_block = spin * TBA_DIM
        DO lat = 0, SUBLATTICES - 2
          current_lat_block = lat * ORBITALS + spin_block + nambu_block
          next_lat_block = current_lat_block + ORBITALS

          !Manually unrolling the double loop over orbitals - the number of orbitals is small and hardcoded anyways.
          Hamiltonian(1 + current_lat_block, 2 + next_lat_block) = Hamiltonian(1 + current_lat_block, 2 + next_lat_block) + Hoppings(hopping_start_idx + 1)
          Hamiltonian(1 + current_lat_block, 3 + next_lat_block) = Hamiltonian(1 + current_lat_block, 3 + next_lat_block) + Hoppings(hopping_start_idx + 2)
          Hamiltonian(2 + current_lat_block, 3 + next_lat_block) = Hamiltonian(2 + current_lat_block, 3 + next_lat_block) + Hoppings(hopping_start_idx + 3)
          Hamiltonian(2 + current_lat_block, 1 + next_lat_block) = Hamiltonian(2 + current_lat_block, 1 + next_lat_block) - Hoppings(hopping_start_idx + 1)
          Hamiltonian(3 + current_lat_block, 1 + next_lat_block) = Hamiltonian(3 + current_lat_block, 1 + next_lat_block) - Hoppings(hopping_start_idx + 2)
          Hamiltonian(3 + current_lat_block, 2 + next_lat_block) = Hamiltonian(3 + current_lat_block, 2 + next_lat_block) - Hoppings(hopping_start_idx + 3)

        END DO
      END DO
    END DO
  END ASSOCIATE

END SUBROUTINE COMPUTE_RASHBA_HOPPING

PURE SUBROUTINE COMPUTE_LAYER_POTENTIAL(Hamiltonian, V_layer, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: V_layer(discretization % SUBLATTICES)
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: nambu, spin, lat, orb, row
  REAL(REAL64) :: sign

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
  REAL(REAL64), INTENT(IN) :: Subband_energies(discretization % SUBBANDS)
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32), INTENT(IN) :: n_band
  INTEGER(INT32) :: nambu, row, i
  REAL(REAL64) :: sign

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
  REAL(REAL64), INTENT(IN) :: E_Fermi
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: nambu, i, row
  REAL(REAL64) :: sign

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
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian of the system that is to be filled
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(discretization % ORBITALS, N_ALL_NEIGHBOURS, SPINS, SPINS, discretization % derived % LAYER_COUPLINGS) !! Superconducting energies
  REAL(REAL64), INTENT(IN) :: kx !! Wavevector in X direction
  REAL(REAL64), INTENT(IN) :: ky !! Wavevector in Y direction
  INTEGER(INT32) :: orb, lat, spin1, spin2, row, col, gamma_lat_index
  INTEGER(INT32) :: spin1_block, spin2_block, spin1_idx, spin2_idx
  INTEGER(INT32) :: current_lat_block, next_lat_block, lat_idx
  INTEGER(INT32) :: neigh

  COMPLEX(REAL64) :: Pairings_nearest(6), Pairings_next(6)

  CALL COMPUTE_NEAREST_PAIRINGS(Pairings_nearest, kx, ky, 6)
  CALL COMPUTE_NEXT_PAIRINGS(Pairings_next, kx, ky, 6)

  !Nearest neighbours pairing
  DO orb = 1, discretization % ORBITALS
    DO spin1 = 0, SPINS - 1
      spin1_block = spin1 * discretization % derived % TBA_DIM
      spin1_idx = spin1 + 1
      DO spin2 = 0, SPINS - 1
        spin2_block = spin2 * discretization % derived % TBA_DIM
        spin2_idx = spin2 + 1
        ! -2, because we have to iterate up to one-before-last sublattice to be able to increment (lat + 1)
        DO lat = 0, discretization % SUBLATTICES - 2
          current_lat_block = lat * discretization % ORBITALS
          next_lat_block = current_lat_block + discretization % ORBITALS

          !Ti1 - Ti2 coupling
          gamma_lat_index = 2 * lat + 2
          row = spin1_block + current_lat_block + orb
          col = discretization % derived % DIM_POSITIVE_K + spin2_block + next_lat_block + orb
          DO neigh = 1, N_NEIGHBOURS
            Hamiltonian(row, col) = Hamiltonian(row, col) + Gamma_SC(orb, neigh, spin1_idx, spin2_idx, gamma_lat_index) * Pairings_nearest(neigh)
          END DO

          !Ti2 - Ti1 coupling
          gamma_lat_index = 2 * lat + 1
          row = spin1_block + next_lat_block + orb
          col = discretization % derived % DIM_POSITIVE_K + spin2_block + current_lat_block + orb
          DO neigh = 1, N_NEIGHBOURS
            Hamiltonian(row, col) = Hamiltonian(row, col) + Gamma_SC(orb, neigh, spin1_idx, spin2_idx, gamma_lat_index) * Pairings_nearest(N_NEIGHBOURS + neigh)
          END DO
        END DO
      END DO
    END DO
  END DO

  !Next nearest neighbours pairing
  DO orb = 1, discretization % ORBITALS
    DO lat = 0, discretization % SUBLATTICES - 1
      current_lat_block = lat * discretization % ORBITALS
      lat_idx = lat + 1
      DO spin1 = 0, SPINS - 1
        spin1_block = spin1 * discretization % derived % TBA_DIM
        spin1_idx = spin1 + 1
        DO spin2 = 0, SPINS - 1
          spin2_block = spin2 * discretization % derived % TBA_DIM
          spin2_idx = spin2 + 1
          row = spin1_block + current_lat_block + orb
          col = discretization % derived % DIM_POSITIVE_K + spin2_block + current_lat_block + orb

          DO neigh = 1, N_NEXT_NEIGHBOURS
            Hamiltonian(row, col) = Hamiltonian(row, col) + Gamma_SC(orb, N_NEIGHBOURS + neigh, spin1_idx, spin2_idx, lat_idx) * Pairings_next(neigh)
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_SC

PURE RECURSIVE SUBROUTINE COMPUTE_HUBBARD(Hamiltonian, Charge_dens, U_HUB, V_HUB, discretization)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: U_HUB, V_HUB
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K)
  INTEGER(INT32) :: orb, lat, orb_prime, spin, nambu, row
  REAL(REAL64) :: sign

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
  REAL(REAL64), INTENT(IN) :: B(3)
  REAL(REAL64), INTENT(IN) :: g_factor
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM)
  INTEGER(INT32) :: i, spin, nambu, row, col
  REAL(REAL64) :: sign_nambu, sign_spin

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

  !H(k)
  DO i = 1, discretization % derived % TBA_DIM
    row = i
    col = discretization % derived % TBA_DIM + i
    !B_x and B_y terms
    Hamiltonian(row, col) = Hamiltonian(row, col) + 0.5d0 * muB * g_factor * (B(1) - imag * B(2))
  END DO

  !-H*(-k)
  DO i = 1, discretization % derived % TBA_DIM
    row = discretization % derived % DIM_POSITIVE_K + i
    col = discretization % derived % DIM_POSITIVE_K + discretization % derived % TBA_DIM + i
    !B_x and B_y terms
    Hamiltonian(row, col) = Hamiltonian(row, col) - 0.5d0 * muB * g_factor * CONJG(B(1) - imag * B(2))
  END DO

END SUBROUTINE COMPUTE_ZEEMAN

PURE RECURSIVE SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING(B, Hamiltonian, discretization)
  !! Computes L \cdot B coupling, taking into account d orbitals
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: B(3) !! Magnetic field
  COMPLEX(REAL64), INTENT(INOUT) :: Hamiltonian(discretization % derived % DIM, discretization % derived % DIM) !! Hamiltonian to be updated

  INTEGER(INT32) :: i_orb, j_orb, lat, spin, nambu, row, col
  REAL(REAL64) :: sign_nambu

  !H(k)
  DO spin = 0, SPINS - 1
    DO lat = 0, discretization % SUBLATTICES - 1
      DO i_orb = 1, discretization % ORBITALS
        DO j_orb = i_orb, discretization % ORBITALS
          row = spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + i_orb
          col = spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + j_orb
          Hamiltonian(row, col) = Hamiltonian(row, col) + muB * &
          & (B(1) * L_x(i_orb, j_orb) + B(2) * L_y(i_orb, j_orb) + B(3) * L_z(i_orb, j_orb))
        END DO
      END DO
    END DO
  END DO

  !-H*(-k)
  DO spin = 0, SPINS - 1
    DO lat = 0, discretization % SUBLATTICES - 1
      DO i_orb = 1, discretization % ORBITALS
        DO j_orb = i_orb, discretization % ORBITALS
          row = discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + i_orb
          col = discretization % derived % DIM_POSITIVE_K + spin * discretization % derived % TBA_DIM + lat * discretization % ORBITALS + j_orb
          Hamiltonian(row, col) = Hamiltonian(row, col) - muB * &
          & CONJG((B(1) * L_x(i_orb, j_orb) + B(2) * L_y(i_orb, j_orb) + B(3) * L_z(i_orb, j_orb)))
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE COMPUTE_ORBITAL_MAGNETIC_COUPLING

END MODULE hamiltonians
