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

MODULE local_integrand
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
USE utilities
USE hamiltonians
USE writers
USE types
IMPLICIT NONE
CONTAINS

SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, Charge_dens, k1, k2,&
                                     & Delta_local, Charge_dens_local, discretization, physical_params)
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  TYPE(physical_params_t), INTENT(IN) :: physical_params
  COMPLEX(REAL64), INTENT(IN) :: Hamiltonian_const(discretization % derived % DIM, &
                                             & discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: k1, k2
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(discretization % ORBITALS, &
                                    & N_ALL_NEIGHBOURS, &
                                    & SPINS, &
                                    & SPINS, &
                                    & discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K)

  COMPLEX(REAL64), INTENT(OUT) :: Delta_local(discretization % ORBITALS, &
                                        & N_ALL_NEIGHBOURS, &
                                        & SPINS, &
                                        & SPINS, &
                                        & discretization % derived % LAYER_COUPLINGS)
  REAL(REAL64), INTENT(OUT) :: Charge_dens_local(discretization % derived % DIM_POSITIVE_K)

  COMPLEX(REAL64) :: Hamiltonian(discretization % derived % DIM, &
                            & discretization % derived % DIM)
  COMPLEX(REAL64) :: U_transformation(discretization % derived % DIM, &
                                & discretization % derived % DIM)
  REAL(REAL64) :: Energies(discretization % derived % DIM)
  REAL(REAL64) :: kx, ky

  !Transform from graphene reciprocal lattice to kx and ky
  kx = k1 * COS(k2)
  ky = k1 * SIN(k2)

  Energies(:) = 0.
  Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
  U_transformation(:, :) = CMPLX(0., 0., KIND=REAL64)
  CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, discretization, physical_params)
  CALL COMPUTE_HUBBARD(Hamiltonian, &
                      & Charge_dens, &
                      & physical_params % subband_params % U_HUB, &
                      & physical_params % subband_params % V_HUB, &
                      & discretization)
  CALL COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, discretization)

  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle

  Hamiltonian(:, :) = 0.5 * (Hamiltonian_const + Hamiltonian)
  !U_transformation(:,:) = Hamiltonian(:,:)
  !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:), DIM)
  ! CALL PRINT_HAMILTONIAN(Hamiltonian(:,:), DIM)
  ! STOP 'Hamiltonian printed'

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian, Energies, U_transformation, discretization % derived % DIM)
  !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U

  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Delta_local = CMPLX(0., 0., KIND=REAL64)
  !Self - consistent delta calculation
  CALL ACCUMULATE_NEAREST_NEIGHBOURS_DELTA(Delta_local, physical_params % subband_params % J_SC_tensor, U_transformation, Energies, kx, ky, discretization, physical_params % external % T)
  CALL ACCUMULATE_NEXT_NEIGHBOURS_DELTA(Delta_local, physical_params % subband_params % J_SC_NNN_tensor, U_transformation, Energies, kx, ky, discretization, physical_params % external % T)

  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Charge_dens_local = 0.
  !Charge density calculation
  CALL ACCUMULATE_CHARGE_DENSITY(Charge_dens_local, U_transformation, Energies, discretization, physical_params % external % T)

  !Multiplication by the Jacobian
  Delta_local = Delta_local * k1
  Charge_dens_local = Charge_dens_local * k1
END SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA

!dir$ attributes forceinline :: ACCUMULATE_NEAREST_NEIGHBOURS_DELTA
PURE SUBROUTINE ACCUMULATE_NEAREST_NEIGHBOURS_DELTA(Delta, J_tensor, U, Energies, kx, ky, discretization, T)
  !! This subroutine computes integrand
  !! <c_{kl\sigma_1} c_{kl\sigma_2}> * exp(i \vec{k} \vec{\delta_{ij}})
  !! For i,j sites being nearest neighbours.
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  COMPLEX(REAL64), INTENT(INOUT) :: Delta(discretization % ORBITALS, &
                                    & N_ALL_NEIGHBOURS, &
                                    & SPINS, &
                                    & SPINS, &
                                    discretization % derived % LAYER_COUPLINGS) !! Accumulator for integrand
  REAL(REAL64), INTENT(IN) :: J_tensor(SPINS, SPINS, SPINS, SPINS)
  COMPLEX(REAL64), INTENT(IN) :: U(discretization % derived % DIM, &
                             & discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64), INTENT(IN) :: Energies(discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), INTENT(IN) :: kx, ky !! Wavevector coordinates
  REAL(REAL64), INTENT(IN) :: T !! Temperature

  INTEGER(INT32) :: orb, lat, spin1, spin2 !! Degrees of freedom of the Hamiltonian
  INTEGER(INT32) :: spin3, spin4 !! Spin degrees of freedom integrated-out in the mean-field approach
  INTEGER(INT32) :: neigh !! neighbor index for phases array
  INTEGER(INT32) :: n !! Index for summation of U transformation matrix
  INTEGER(INT32) :: row, col !! Postion in the Hamiltonian based on degrees of freedom indeces
  INTEGER(INT32) :: row_inv, col_inv !! Postion in the Hamiltonian for opposite lattice hopping based on degrees of freedom indeces
  INTEGER(INT32) :: lat_block_1, lat_block_2 !! Shift in hamiltonian indeces from the lattice degree of freedom
  INTEGER(INT32) :: spin3_block !! Shift in hamiltonian indeces from the spin degree of freedom
  INTEGER(INT32) :: spin4_block !! Shift in hamiltonian indeces from the spin degree of freedom
  INTEGER(INT32) :: lat_idx_inv
  REAL(REAL64) :: occupation_electron, occupation_hole !! Occupation of a given (n-th) energy level
  REAL(REAL64) :: j_elem !! Energy of nearest neighbour superconducting pairing
  COMPLEX(REAL64) :: average_pairing, average_energy
  COMPLEX(REAL64) :: average_pairing_inv, average_energy_inv
  COMPLEX(REAL64) :: phases(3) !! Phase factors for nearest neighbours

  phases(1) = pairing_1(ky)
  phases(2) = pairing_2(kx, ky)
  phases(3) = pairing_3(kx, ky)

  !TODO: Think about this spin1,2,3,4 summation and which index has to be taken into account in row/col and which in Delta() indexing.
  DO orb = 1, discretization % ORBITALS
    !Electrons
    DO n = 1, discretization % derived % DIM_POSITIVE_K
      occupation_electron = fd_distribution(Energies(n), 0d0, T)
      occupation_hole = 1.0 - fd_distribution(-Energies(discretization % derived % DIM_POSITIVE_K + n), 0d0, T)
      DO lat = 1, discretization % derived % LAYER_COUPLINGS, 2
        lat_block_1 = (lat / 2) * discretization % ORBITALS
        lat_block_2 = ((lat + 1) / 2) * discretization % ORBITALS
        lat_idx_inv = lat + 1
        DO spin3 = 1, SPINS
          spin3_block = (spin3 - 1) * discretization % derived % TBA_DIM
          row = orb + spin3_block + lat_block_1 + discretization % derived % DIM_POSITIVE_K
          row_inv = orb + spin3_block + lat_block_2 + discretization % derived % DIM_POSITIVE_K
          DO spin4 = 1, SPINS
            spin4_block = (spin4 - 1) * discretization % derived % TBA_DIM
            col = orb + spin4_block + lat_block_2
            col_inv = orb + spin4_block + lat_block_1

            average_pairing = CONJG(U(row, n)) * U(col, n) * occupation_electron + &
            & CONJG(U(row, discretization % derived % DIM_POSITIVE_K + n)) * U(col, discretization % derived % DIM_POSITIVE_K + n) * occupation_hole

            average_pairing_inv = CONJG(U(row_inv, n)) * U(col_inv, n) * occupation_electron + &
            & CONJG(U(row_inv, discretization % derived % DIM_POSITIVE_K + n)) * U(col_inv, discretization % derived % DIM_POSITIVE_K + n) * occupation_hole

            DO spin1 = 1, SPINS
              DO spin2 = 1, SPINS
                j_elem = J_tensor(spin1, spin2, spin3, spin4)
                average_energy = j_elem * average_pairing
                average_energy_inv = j_elem * average_pairing_inv
                DO neigh = 1, N_NEIGHBOURS
                  Delta(orb, neigh, spin1, spin2, lat) = Delta(orb, neigh, spin1, spin2, lat) + average_energy * phases(neigh)
                  Delta(orb, neigh, spin1, spin2, lat_idx_inv) = Delta(orb, neigh, spin1, spin2, lat_idx_inv) + average_energy_inv * CONJG(phases(neigh))
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ACCUMULATE_NEAREST_NEIGHBOURS_DELTA

!dir$ attributes forceinline :: ACCUMULATE_NEXT_NEIGHBOURS_DELTA
PURE SUBROUTINE ACCUMULATE_NEXT_NEIGHBOURS_DELTA(Delta, J_tensor, U, Energies, kx, ky, discretization, T)
  !! This subroutine computes integrand
  !! <c_{kl\sigma_1} c_{kl\sigma_2}> * exp(i \vec{k} \vec{\delta_{ij}})
  !! For i,j sites being next-nearest neighbours.
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  COMPLEX(REAL64), INTENT(INOUT) :: Delta(discretization % ORBITALS, &
                                    & N_ALL_NEIGHBOURS, &
                                    & SPINS, &
                                    & SPINS, &
                                    & discretization % derived % LAYER_COUPLINGS) !! Accumulator for integrand
  REAL(REAL64), INTENT(IN) :: J_tensor(SPINS, SPINS, SPINS, SPINS)
  COMPLEX(REAL64), INTENT(IN) :: U(discretization % derived % DIM, &
                             & discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64), INTENT(IN) :: Energies(discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), INTENT(IN) :: kx, ky !! Wavevector coordinates
  REAL(REAL64), INTENT(IN) :: T !! Temperature

  INTEGER(INT32) :: orb, lat, spin1, spin2 !! Degrees of freedom of the Hamiltonian
  INTEGER(INT32) :: spin3, spin4 !! Spin degrees of freedom integrated-out in the mean-field approach
  INTEGER(INT32) :: neigh !! neighbor index for phases array
  INTEGER(INT32) :: n !! Index for summation of U transformation matrix
  INTEGER(INT32) :: row, col !! Postion in the Hamiltonian based on degrees of freedom indeces
  INTEGER(INT32) :: lat_block !! Shift in hamiltonian indeces from the lattice degree of freedom
  INTEGER(INT32) :: spin3_block !! Shift in hamiltonian indeces from the spin degree of freedom
  INTEGER(INT32) :: spin4_block !! Shift in hamiltonian indeces from the spin degree of freedom
  INTEGER(INT32) :: lat_idx
  REAL(REAL64) :: occupation_electron, occupation_hole !! Occupation of a given (n-th) energy level
  REAL(REAL64) :: j_elem !! Energy of nearest neighbour superconducting pairing
  COMPLEX(REAL64) :: average_pairing, average_energy
  COMPLEX(REAL64) :: phases(6) !! Phase factors for next nearest neighbours

  phases(1) = CONJG(pairing_nnn_1(kx))
  phases(2) = CONJG(pairing_nnn_2(kx, ky))
  phases(3) = CONJG(pairing_nnn_3(kx, ky))
  phases(4) = CONJG(pairing_nnn_4(kx))
  phases(5) = CONJG(pairing_nnn_5(kx, ky))
  phases(6) = CONJG(pairing_nnn_6(kx, ky))

  DO orb = 1, discretization % ORBITALS
    DO n = 1, discretization % derived % DIM_POSITIVE_K
      occupation_electron = fd_distribution(Energies(n), 0d0, T)
      occupation_hole = 1.0 - fd_distribution(-Energies(discretization % derived % DIM_POSITIVE_K + n), 0d0, T)
      !Up - down Ti1 - Ti1 delta, Ti2 - Ti2 delta
      !No conjugation in phase factor, since next nearest neighbours have the same relative positions in both sublattices
      DO lat = 0, discretization % SUBLATTICES - 1
        lat_block = lat * discretization % ORBITALS
        lat_idx = lat + 1
        DO spin3 = 1, SPINS
          row = orb + lat_block + (spin3 - 1) * discretization % derived % TBA_DIM + discretization % derived % DIM_POSITIVE_K
          DO spin4 = 1, SPINS

            col = orb + lat_block + (spin4 - 1) * discretization % derived % TBA_DIM

            average_pairing = CONJG(U(row, n)) * U(col, n) * occupation_electron + &
            & CONJG(U(row, discretization % derived % DIM_POSITIVE_K + n)) * U(col, discretization % derived % DIM_POSITIVE_K + n) * occupation_hole

            DO spin1 = 1, SPINS
              DO spin2 = 1, SPINS

                ! Minus because we need the potential to be attractive. This way we can pass positive values in input.nml
                j_elem = J_tensor(spin1, spin2, spin3, spin4)
                average_energy = j_elem * average_pairing

                DO neigh = 1, N_ALL_NEIGHBOURS - N_NEIGHBOURS
                  Delta(orb, N_NEIGHBOURS + neigh, spin1, spin2, lat_idx) = Delta(orb, N_NEIGHBOURS + neigh, spin1, spin2, lat_idx) + &
                  & average_energy * phases(neigh)
                END DO

              END DO
            END DO

          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE ACCUMULATE_NEXT_NEIGHBOURS_DELTA

SUBROUTINE ACCUMULATE_CHARGE_DENSITY(Charge, U, Energies, discretization, T)
  !! This subroutine computes integrand
  !! <c_{kl\sigma}^\dag c_{kl\sigma}>
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(INOUT) :: Charge(discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64), INTENT(IN) :: U(discretization % derived % DIM, discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64), INTENT(IN) :: Energies(discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), INTENT(IN) :: T !! Temperature

  INTEGER(INT32) :: n, m !! Index for summation of U transformation matrix

  DO m = 1, discretization % derived % DIM_POSITIVE_K
    DO n = 1, discretization % derived % DIM_POSITIVE_K
      Charge(m) = Charge(m) + REAL(U(m, n) * CONJG(U(m, n)), KIND=8) * fd_distribution(Energies(n), 0d0, T) + &
      & REAL(U(m, discretization % derived % DIM_POSITIVE_K + n) * CONJG(U(m, discretization % derived % DIM_POSITIVE_K + n)), KIND=8) * &
      & (1.-fd_distribution(-Energies(discretization % derived % DIM_POSITIVE_K + n), 0d0, T))
    END DO
  END DO

END SUBROUTINE ACCUMULATE_CHARGE_DENSITY

END MODULE local_integrand
