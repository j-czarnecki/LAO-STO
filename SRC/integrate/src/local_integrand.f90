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
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                        & discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64), INTENT(OUT) :: Delta_local(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                            & discretization % derived % DIM_POSITIVE_K, &
                                            & discretization % derived % DIM_POSITIVE_K)
#else
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64), INTENT(OUT) :: Delta_local(discretization % derived % DIM_POSITIVE_K, &
                                            & discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64) :: Hamiltonian_electron(discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64) :: Hamiltonian_hole(discretization % derived % DIM_POSITIVE_K, &
                                    & discretization % derived % DIM_POSITIVE_K)
  COMPLEX(REAL64) :: Hamiltonian_hole_reversed(discretization % derived % DIM_POSITIVE_K, &
                                             & discretization % derived % DIM_POSITIVE_K)
  REAL(REAL64) :: Energies_electron(discretization % derived % DIM_POSITIVE_K)
  REAL(REAL64) :: Energies_hole(discretization % derived % DIM_POSITIVE_K)
  INTEGER(INT32) :: i, j
#endif
  COMPLEX(REAL64) :: U_transformation(discretization % derived % DIM, &
                                & discretization % derived % DIM)
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K)
  REAL(REAL64), INTENT(OUT) :: Charge_dens_local(discretization % derived % DIM_POSITIVE_K)

  COMPLEX(REAL64) :: Hamiltonian(discretization % derived % DIM, &
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

#ifndef BAND_BASIS
  CALL COMPUTE_HUBBARD(Hamiltonian, &
                      & Charge_dens, &
                      & physical_params % subband_params % U_HUB, &
                      & physical_params % subband_params % V_HUB, &
                      & discretization)
  CALL COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle

  Hamiltonian(:, :) = 0.5 * (Hamiltonian_const + Hamiltonian)

#else
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle
  Hamiltonian = 0.5 * (Hamiltonian_const + Hamiltonian)
  Hamiltonian_electron = Hamiltonian(:discretization % derived % DIM_POSITIVE_K, &
                                   & :discretization % derived % DIM_POSITIVE_K)
  Hamiltonian_hole = Hamiltonian(discretization % derived % DIM_POSITIVE_K + 1:, &
                               & discretization % derived % DIM_POSITIVE_K + 1:)
  CALL DIAGONALIZE_HERMITIAN(Hamiltonian_electron, Energies_electron, discretization % derived % DIM_POSITIVE_K)
  CALL DIAGONALIZE_HERMITIAN(Hamiltonian_hole, Energies_hole, discretization % derived % DIM_POSITIVE_K)

  Hamiltonian = CMPLX(0., 0., KIND=REAL64)
  CALL COMPUTE_HUBBARD(Hamiltonian, &
                      & Charge_dens, &
                      & physical_params % subband_params % U_HUB, &
                      & physical_params % subband_params % V_HUB, &
                      & discretization)
  CALL COMPUTE_SC_BAND(Hamiltonian, kx, ky, Gamma_SC, discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle
  Hamiltonian = 0.5 * Hamiltonian
  DO i = 1, discretization % derived % DIM_POSITIVE_K
    Hamiltonian(i, i) = Hamiltonian(i, i) + Energies_electron(i)
    ! Since ZHEEV sorts energies in ascending order, we have to fill hole energies in reverse order.
    ! This guarantees that the SC pairing block is diagonal (not antidiagonal).
    Hamiltonian(discretization % derived % DIM_POSITIVE_K + i, discretization % derived % DIM_POSITIVE_K + i) = &
      & Hamiltonian(discretization % derived % DIM_POSITIVE_K + i, discretization % derived % DIM_POSITIVE_K + i) + &
      & Energies_hole(discretization % derived % DIM_POSITIVE_K + 1 - i)
  END DO
#endif
  CALL DIAGONALIZE_GENERALIZED(Hamiltonian, Energies, U_transformation, discretization % derived % DIM)
  !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U

  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Delta_local = CMPLX(0., 0., KIND=REAL64)
  !Self - consistent delta calculation
#ifndef BAND_BASIS
  CALL ACCUMULATE_DELTA_REAL_SPACE(Delta_local, physical_params % subband_params % J_tensor % Values, &
    & physical_params % subband_params % J_tensor % Column_indices, &
    & physical_params % subband_params % J_tensor % Row_indices, &
    & U_transformation, Energies, kx, ky, discretization, &
    & physical_params % subband_params % J_tensor % n_nonzero, physical_params % external % T)
#else
  CALL ACCUMULATE_DELTA_K_SPACE(Delta_local, physical_params % subband_params % J_tensor % Values, &
    & physical_params % subband_params % J_tensor % Column_indices, &
    & physical_params % subband_params % J_tensor % Row_indices, &
    & U_transformation, Energies, kx, ky, discretization, &
    & physical_params % subband_params % J_tensor % n_nonzero, physical_params % external % T)
#endif
  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Charge_dens_local = 0.
  !Charge density calculation
  CALL ACCUMULATE_CHARGE_DENSITY(Charge_dens_local, U_transformation, Energies, discretization, physical_params % external % T)
  !Multiplication by the Jacobian
  Delta_local = Delta_local * k1
  Charge_dens_local = Charge_dens_local * k1
END SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA

!dir$ attributes forceinline :: ACCUMULATE_DELTA_REAL_SPACE
SUBROUTINE ACCUMULATE_DELTA_REAL_SPACE(Delta, Values_j_tensor, Column_indeces_j_tensor, Row_indeces_j_tensor,&
                                             & U, Energies, kx, ky, discretization, nonzero_j_tensor, T)
  !! This subroutine computes integrand
  !! <c_{kl\sigma_1} c_{kl\sigma_2}> * exp(i \vec{k} \vec{\delta_{ij}})
  !! For i,j sites being nearest neighbours.
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization
  INTEGER(INT32), INTENT(IN) :: nonzero_j_tensor
  COMPLEX(REAL64), INTENT(INOUT) :: Delta(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                    & discretization % derived % DIM_POSITIVE_K, &
                                    & discretization % derived % DIM_POSITIVE_K) !! Accumulator for integrand
  REAL(REAL64), INTENT(IN) :: Values_j_tensor(nonzero_j_tensor) !! Nonzero values of the interaction tensor
  INTEGER(INT32), INTENT(IN) :: Column_indeces_j_tensor(nonzero_j_tensor) !! Column indeces of nonzero elements of the interaction tensor
  INTEGER(INT32), INTENT(IN) :: Row_indeces_j_tensor(discretization % derived % DIM_POSITIVE_K**2 + 1) !! Indeces of elements in
                                                                                                       !! Values_j_tensor/Column_indeces_j_tensor that start
                                                                                                       !! a new row.
  COMPLEX(REAL64), INTENT(IN) :: U(discretization % derived % DIM, &
                                 & discretization % derived % DIM) !! Unitary matrix that diagonalizes the Hamiltonian - from ZGEEV
  REAL(REAL64), INTENT(IN) :: Energies(discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), INTENT(IN) :: kx, ky !! Wavevector coordinates
  REAL(REAL64), INTENT(IN) :: T !! Temperature

  INTEGER(INT32) :: neigh !! neighbor index for phases array
  INTEGER(INT32) :: n !! Index for summation of U transformation matrix
  INTEGER(INT32) :: row, col !! Postion in the Hamiltonian based on degrees of freedom indeces
  REAL(REAL64) :: occupation_electron, occupation_hole !! Occupation of a given (n-th) energy level
  COMPLEX(REAL64) :: average_pairing, average_energy
  COMPLEX(REAL64) :: Phases(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS) !! Phase factors for subsequent neighbours

  INTEGER(INT32) :: i, j
  INTEGER(INT32) :: Dematricized_indeces_first(2) !! Two-tuple of band indeces (\alpha, \beta) from the row index of interaction tensor.
  INTEGER(INT32) :: Dematricized_indeces_second(2) !! Two-tuple of band indeces (\gamma, \delta) from the column index of interaction tensor.

  CALL COMPUTE_NEAREST_PAIRINGS(Phases(:N_NEAREST_NEIGHBOURS), kx, ky, N_NEAREST_NEIGHBOURS)
  CALL COMPUTE_NEXT_PAIRINGS(Phases(N_NEAREST_NEIGHBOURS + 1:), kx, ky, N_NEXT_NEIGHBOURS)

  DO n = 1, discretization % derived % DIM_POSITIVE_K
    occupation_electron = fd_distribution(Energies(n), 0d0, T)
    occupation_hole = 1.0 - fd_distribution(-Energies(discretization % derived % DIM_POSITIVE_K + n), 0d0, T)

    DO i = 1, discretization % derived % DIM_POSITIVE_K**2 !! Loop over all rows of tensor
      Dematricized_indeces_first = get_dematricized_indeces(i, discretization % derived % DIM_POSITIVE_K)
      DO j = Row_indeces_j_tensor(i), Row_indeces_j_tensor(i + 1) - 1 !! Loop over all columns of tensor
        Dematricized_indeces_second = get_dematricized_indeces(Column_indeces_j_tensor(j), discretization % derived % DIM_POSITIVE_K)
        row = Dematricized_indeces_second(1) + discretization % derived % DIM_POSITIVE_K
        col = Dematricized_indeces_second(2)
        average_pairing = CONJG(U(row, n)) * U(col, n) * occupation_electron + &
          & CONJG(U(row, discretization % derived % DIM_POSITIVE_K + n)) * U(col, discretization % derived % DIM_POSITIVE_K + n) * occupation_hole
        average_energy = Values_j_tensor(j) * average_pairing
        !! Calculating for each nearest neighbour
        DO neigh = 1, N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS
          Delta(neigh, Dematricized_indeces_first(1), Dematricized_indeces_first(2)) = &
            & Delta(neigh, Dematricized_indeces_first(1), Dematricized_indeces_first(2)) + &
            & average_energy * Phases(neigh)
        END DO
      END DO
    END DO ! Loop over elements of pairing matrix
  END DO ! Transformation matrix loop

END SUBROUTINE ACCUMULATE_DELTA_REAL_SPACE

!SUBROUTINES for interaction tensor expressed in band basis
!dir$ attributes forceinline :: ACCUMULATE_DELTA_K_SPACE
SUBROUTINE ACCUMULATE_DELTA_K_SPACE(Delta, Values_j_tensor, Column_indeces_j_tensor, Row_indeces_j_tensor, &
                                     & U, Energies, kx, ky, discretization, nonzero_j_tensor, T)
  !! This subroutine computes integrand <c_{-k \gamma}^\dag c_{k \delta}>,
  !! where \gamma and \delta are band indeces.
  !! It is calculated based on assumption that interaction tensor V^{\alpha \beta \gamma \delta}
  !! is independent of the wavevector.
  !! As a result it updates the accumulator for anomalous averages (Delta).
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization !! Discretization parameters
  INTEGER(INT32), INTENT(IN) :: nonzero_j_tensor !! Number of non-zero elements in the interaction tensor
  COMPLEX(REAL64), INTENT(INOUT) :: Delta(discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K) !! Accumulator for integrand
  REAL(REAL64), INTENT(IN) :: Values_j_tensor(nonzero_j_tensor) !! Nonzero values of the interaction tensor
  INTEGER(INT32), INTENT(IN) :: Column_indeces_j_tensor(nonzero_j_tensor) !! Column indeces of nonzero elements of the interaction tensor
  INTEGER(INT32), INTENT(IN) :: Row_indeces_j_tensor(discretization % derived % DIM_POSITIVE_K**2 + 1) !! Indeces of elements in
                                                                                                       !! Values_j_tensor/Column_indeces_j_tensor that start
                                                                                                       !! a new row.
  COMPLEX(REAL64), INTENT(IN) :: U(discretization % derived % DIM, &
                                 & discretization % derived % DIM) !! Unitary matrix that diagonalizes the interacting Hamiltonian
  REAL(REAL64), INTENT(IN) :: Energies(discretization % derived % DIM) !! Energies for given wavevector
  REAL(REAL64), INTENT(IN) :: kx, ky !! Wavevector coordinates
  REAL(REAL64), INTENT(IN) :: T !! Temperature

  INTEGER(INT32) :: i, j, k, l !! Band indeces
  INTEGER(INT32) :: Dematricized_indeces_first(2) !! Two-tuple of band indeces (\alpha, \beta) from the row index of interaction tensor.
  INTEGER(INT32) :: Dematricized_indeces_second(2) !! Two-tuple of band indeces (\gamma, \delta) from the column index of interaction tensor.
  INTEGER(INT32) :: n !! Column from unitary matrix that diagonalizes full Hamiltonian with interaction.
  INTEGER(INT32) :: row, col !! Postion in the Hamiltonian/unitary matrix.
  REAL(REAL64) :: occupation_electron, occupation_hole !! Occupation of a given (n-th) energy level.
  COMPLEX(REAL64) :: average_pairing !! Average pairing for a given (n-th) energy level.
  COMPLEX(REAL64) :: average_energy !! Average pairing multiplied by tensor elemment i.e. "energy".

  DO n = 1, discretization % derived % DIM_POSITIVE_K
    occupation_electron = fd_distribution(Energies(n), 0d0, T)
    occupation_hole = 1.0 - fd_distribution(-Energies(discretization % derived % DIM_POSITIVE_K + n), 0d0, T)

    DO i = 1, discretization % derived % DIM_POSITIVE_K**2 !! Loop over all rows of tensor
      Dematricized_indeces_first = get_dematricized_indeces(i, discretization % derived % DIM_POSITIVE_K)
      DO j = Row_indeces_j_tensor(i), Row_indeces_j_tensor(i + 1) - 1 !! Loop over all columns of tensor
        Dematricized_indeces_second = get_dematricized_indeces(Column_indeces_j_tensor(j), discretization % derived % DIM_POSITIVE_K)
        row = Dematricized_indeces_second(1) + discretization % derived % DIM_POSITIVE_K
        col = Dematricized_indeces_second(2)
        average_pairing = CONJG(U(row, n)) * U(col, n) * occupation_electron + &
          & CONJG(U(row, discretization % derived % DIM_POSITIVE_K + n)) * U(col, discretization % derived % DIM_POSITIVE_K + n) * occupation_hole
        average_energy = Values_j_tensor(j) * average_pairing

        Delta(Dematricized_indeces_first(1), Dematricized_indeces_first(2)) = &
          & Delta(Dematricized_indeces_first(1), Dematricized_indeces_first(2)) + average_energy
      END DO
    END DO ! Loop over elements of pairing matrix
  END DO ! Transformation matrix loop

END SUBROUTINE ACCUMULATE_DELTA_K_SPACE

!dir$ attributes forceinline :: ACCUMULATE_CHARGE_DENSITY
SUBROUTINE ACCUMULATE_CHARGE_DENSITY(Charge, U, Energies, discretization, T)
  !! This subroutine computes integrand
  !! <c_{kl\sigma}^\dag c_{kl\sigma}>
  IMPLICIT NONE
  TYPE(discretization_t), INTENT(IN) :: discretization !! Discretizatin parameters
  REAL(REAL64), INTENT(INOUT) :: Charge(discretization % derived % DIM_POSITIVE_K) !! Accumulator for integrand of charge density
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
