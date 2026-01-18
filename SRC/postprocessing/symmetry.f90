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

#include "macros_def.f90"
MODULE symmetry
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE hamiltonians
USE parameters
USE utilities
USE writers
USE reader
USE local_integrand
USE self_consistency
USE logger
USE types
IMPLICIT NONE
CONTAINS

SUBROUTINE CALCULATE_GAMMA_K(gamma)
  TYPE(post_gamma_k_t), INTENT(IN) :: gamma
  TYPE(sc_input_params_t) :: sc_input
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian_const(:, :) !! k-indepndent and band-independent part of the Hamiltonian
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian_const_band(:, :) !! k-independent, band-dependent Hamiltonian
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian_dummy(:, :) !! This is only used to get matrix elements gamma that show up in the Hamiltonian
#ifndef BAND_BASIS
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :) !! Supercondicting pairings read from a simulation
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_K(:, :, :, :) !! Superconducting pairing, determined at given k point
  COMPLEX(REAL64), ALLOCATABLE :: Delta_local(:, :, :, :) !! Delta (pairing amplitudes) for given k point, integrand
#else
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :)
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_K(:, :, :) !! Superconducting pairing, determined at given k point
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_K_orig_basis(:, :) !! Superconducting pairing, determined at given k point
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian(:, :) !! Full Hamiltonian at given k point
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian_electron(:, :), Hamiltonian_hole(:, :)
  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian_hole_reversed(:, :)
  REAL(REAL64), ALLOCATABLE :: Energies_electron(:), Energies_hole(:)
  COMPLEX(REAL64), ALLOCATABLE :: Delta_local(:, :, :) !! Delta (pairing amplitudes) for given k point, integrand
#endif
  REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :) !! Charge density read from a simulation
  REAL(REAL64), ALLOCATABLE :: Charge_dens_local(:, :) !! Charge density for given k point, integrand

  REAL(REAL64) :: k1, k2, kx, ky, kx_minus, ky_minus, dkx, dky
  INTEGER(INT32) :: i, j, n, m, band
  INTEGER(INT32) :: kx_steps, ky_steps
  INTEGER(INT32) :: orb, neigh, spin1, spin2, layer, lat, file_count
  INTEGER(INT32) :: orb_prime, band_prime
  INTEGER(INT32) :: row, col, row_inverse, col_inverse, row_nnn, col_nnn
  INTEGER(INT32) :: gamma_lat_index, gamma_spin_index
  INTEGER(INT32) :: i_band, j_band
  INTEGER(INT32) :: idx
  COMPLEX(REAL64) :: phase

  COMPLEX(REAL64) :: gamma_nn_12, gamma_nn_21, gamma_nnn

  CHARACTER(LEN=200) :: filename
  COMPLEX(REAL64), ALLOCATABLE :: File_unit_mapping(:, :, :)

  LOGICAL :: file_exists

  REAL(REAL64) :: Brillouin_zone_vertices(6, 2)

  Brillouin_zone_vertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  Brillouin_zone_vertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  dkx = 3 * KX_MAX / gamma % Nk_points
  dky = 3 * KY_MAX / gamma % Nk_points

  kx_steps = INT(gamma % Nk_points)
  ky_steps = INT(gamma % Nk_points)

  CALL GET_INPUT(TRIM(gamma % path)//"input.nml", sc_input)
  ASSOCIATE (SUBLATTICES => sc_input % discretization % SUBLATTICES, &
         & SUBBANDS => sc_input % discretization % SUBBANDS, &
         & ORBITALS => sc_input % discretization % ORBITALS, &
         & TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
         & DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
         & DIM => sc_input % discretization % derived % DIM, &
         & LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)

    ALLOCATE (Hamiltonian_const(DIM, DIM))
    ALLOCATE (Hamiltonian_const_band(DIM, DIM))
    ALLOCATE (Hamiltonian_dummy(DIM, DIM))
#ifndef BAND_BASIS
    ALLOCATE (Gamma_SC(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Gamma_K(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Delta_local(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#else
    ALLOCATE (Gamma_SC(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Gamma_K(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Gamma_K_orig_basis(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (Delta_local(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Hamiltonian_electron(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (Hamiltonian_hole(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (Hamiltonian_hole_reversed(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (Hamiltonian(DIM, DIM))
    ALLOCATE (Energies_electron(DIM_POSITIVE_K))
    ALLOCATE (Energies_hole(DIM_POSITIVE_K))
#endif
    ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Charge_dens_local(DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (File_unit_mapping(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
  END ASSOCIATE

  Hamiltonian_const = CMPLX(0., 0., KIND=REAL64)
  Hamiltonian_const_band = CMPLX(0., 0., KIND=REAL64)
  Gamma_SC = CMPLX(0., 0., KIND=REAL64)
  Charge_dens = 0.

  CALL GET_SAFE_GAMMA_SC(Gamma_SC, gamma % path, sc_input % discretization)
  CALL GET_SAFE_CHARGE_DENS(Charge_dens, gamma % path, sc_input % discretization)

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

  !Opening all files I will write gammas to and create a mapping of file units to appropriate names
  file_count = 10
  DO i_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
    DO j_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
      DO band = 1, sc_input % discretization % SUBBANDS
        File_unit_mapping(i_band, j_band, band) = file_count
        WRITE (filename, "(3(A, I0))") "GammaK_i", i_band, "_j", j_band, "_band", band
        OPEN (unit=file_count, FILE=TRIM(gamma % path)//"OutputData/"//TRIM(filename)//".dat", FORM="FORMATTED", ACTION="WRITE")
        WRITE (file_count, '(A)') "#kx[1/a]   ky[1/a]    Re(H_{i,j})[meV]     Im(H_{i,j})[meV]    Re(Delta_neigh1)[meV]   Im(Delta_neigh1)[meV]   Re(Delta_neigh2)[meV] ..."
        file_count = file_count + 1
      END DO
    END DO
  END DO

  !$omp parallel private(kx, ky, kx_minus, ky_minus, k1, k2, band, Hamiltonian_const_band, Charge_dens_local, Delta_local,&
  !$omp&                 Gamma_K, orb, spin1, spin2, n, neigh, lat, orb_prime, band_prime, layer, file_count, row, col, &
  !$omp&                 gamma_spin_index, gamma_lat_index, Hamiltonian_dummy, gamma_nn_12, gamma_nn_21, gamma_nnn, idx, phase &
#ifdef BAND_BASIS
  !$omp&             ,   Hamiltonian, Hamiltonian_electron, Hamiltonian_hole, Energies_electron, Energies_hole, &
  !$omp&                 Hamiltonian_hole_reversed, Gamma_K_orig_basis &
#endif
  !$omp&             )

  !$omp do
  DO i = -kx_steps, kx_steps
    DO j = -ky_steps, ky_steps
      kx = i * dkx
      ky = j * dky

      DO band = 1, sc_input % discretization % SUBBANDS
        Hamiltonian_const_band = Hamiltonian_const
        CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band, sc_input % physical % subband_params % Subband_energies, sc_input % discretization)

        k1 = SQRT(kx**2 + ky**2)
        k2 = ATAN2(ky, kx)
#ifndef BAND_BASIS
        CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const_band, &
                                        & Gamma_SC(:, :, :, band), &
                                        & Charge_dens(:, band), &
                                        & k1, &
                                        & k2, &
                                        & Delta_local(:, :, :, band), &
                                        & Charge_dens_local(:, band), &
                                        & sc_input % discretization, &
                                        & sc_input % physical)
#else
        CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const_band, &
                                        & Gamma_SC(:, :, band), &
                                        & Charge_dens(:, band), &
                                        & k1, &
                                        & k2, &
                                        & Delta_local(:, :, band), &
                                        & Charge_dens_local(:, band), &
                                        & sc_input % discretization, &
                                        & sc_input % physical)
#endif
      END DO

      !Write result for given k point to file
      DO band = 1, sc_input % discretization % SUBBANDS
        Hamiltonian_dummy = CMPLX(0.0d0, 0.0d0, KIND=REAL64)
#ifndef BAND_BASIS
        CALL COMPUTE_SC(Hamiltonian_dummy, kx, ky, Gamma_SC(:, :, :, band), sc_input % discretization)
        DO i_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          DO j_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
            file_count = File_unit_mapping(i_band, j_band, band)
            WRITE (file_count, '(2E15.5, *(2E15.5))') kx, ky, &
              & REAL(Hamiltonian_dummy(i_band, sc_input % discretization % derived % DIM_POSITIVE_K + j_band)) / meV2au, &
              & AIMAG(Hamiltonian_dummy(i_band, sc_input % discretization % derived % DIM_POSITIVE_K + j_band)) / meV2au, &
              & (REAL(Delta_local(neigh, i_band, j_band, band)) / meV2au, &
              & AIMAG(Delta_local(neigh, i_band, j_band, band)) / meV2au, &
              & neigh=1, &
              & N_NEXT_NEIGHBOURS + N_NEAREST_NEIGHBOURS)
          END DO
        END DO
#else
        Hamiltonian = Hamiltonian_const_band
        CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
        Hamiltonian = 0.5 * Hamiltonian

        Hamiltonian_electron = Hamiltonian(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                                         & :sc_input % discretization % derived % DIM_POSITIVE_K)
        CALL DIAGONALIZE_HERMITIAN(Hamiltonian_electron, Energies_electron, sc_input % discretization % derived % DIM_POSITIVE_K)
        !! Fix the gauge
        DO n = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          idx = MAXLOC(ABS(Hamiltonian_electron(:, n)), 1)
          phase = Hamiltonian_electron(idx, n) / ABS(Hamiltonian_electron(idx, n))
          Hamiltonian_electron(:, n) = Hamiltonian_electron(:, n) / phase
        END DO

        !! Electron hamiltonian at -k
        kx_minus = -kx
        ky_minus = -ky
        Hamiltonian = Hamiltonian_const_band
        CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx_minus, ky_minus, sc_input % discretization, sc_input % physical)
        Hamiltonian = 0.5 * Hamiltonian

        Hamiltonian_hole = Hamiltonian(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                                     & :sc_input % discretization % derived % DIM_POSITIVE_K)
        CALL DIAGONALIZE_HERMITIAN(Hamiltonian_hole, Energies_hole, sc_input % discretization % derived % DIM_POSITIVE_K)
        !! Fix the gauge
        DO n = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          idx = MAXLOC(ABS(Hamiltonian_hole(:, n)), 1)
          phase = Hamiltonian_hole(idx, n) / ABS(Hamiltonian_hole(idx, n))
          Hamiltonian_hole(:, n) = Hamiltonian_hole(:, n) / phase
        END DO

        CALL COMPUTE_SC_BAND(Hamiltonian_dummy, kx, ky, Gamma_SC(:, :, band), sc_input % discretization)
        ! Transform back to spin-orbital-sublattice basis to get Hamiltonian elements
        Gamma_K(:, :, band) = Hamiltonian_dummy(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                                        & sc_input % discretization % derived % DIM_POSITIVE_K + 1:)
        Gamma_K_orig_basis(:, :) = MATMUL(Hamiltonian_electron, MATMUL((Gamma_K(:, :, band)), TRANSPOSE(Hamiltonian_hole)))

        DO i_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          DO j_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
            file_count = File_unit_mapping(i_band, j_band, band)
            WRITE (file_count, '(2E15.5, *(2E15.5))') kx, ky, &
              & REAL(Gamma_K_orig_basis(i_band, j_band)) / meV2au, &
              & AIMAG(Gamma_K_orig_basis(i_band, j_band)) / meV2au, &
              & (REAL(Delta_local(i_band, j_band, band)) / meV2au, &
              & AIMAG(Delta_local(i_band, j_band, band)) / meV2au, &
              & neigh=1, &
              & N_NEXT_NEIGHBOURS + N_NEAREST_NEIGHBOURS) !! Repeating values to unify format between bases.
          END DO
        END DO

#endif

      END DO
    END DO
  END DO
  !$omp end do
  !$omp end parallel

  !Close all files
  DO i_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
    DO j_band = 1, sc_input % discretization % derived % DIM_POSITIVE_K
      DO band = 1, sc_input % discretization % SUBBANDS
        file_count = File_unit_mapping(i_band, j_band, band)
        CLOSE (file_count)
      END DO
    END DO
  END DO

  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (Hamiltonian_const_band)
  DEALLOCATE (Hamiltonian_dummy)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Gamma_K)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (Delta_local)
  DEALLOCATE (Charge_dens_local)
  DEALLOCATE (File_unit_mapping)
END SUBROUTINE CALCULATE_GAMMA_K

!! DEPRECATED
SUBROUTINE CALCULATE_PROJECTIONS(projections)
  TYPE(post_projections_t), INTENT(IN) :: projections
  TYPE(sc_input_params_t) :: sc_input
  !TODO: Enable spin-triplet projections
  ! CHARACTER(LEN=*), INTENT(IN) :: input_path !! This should be a path to folder where input.nml resides
  ! INTEGER, INTENT(IN) :: n_r_points !! Number of steps in radial direction
  ! INTEGER, INTENT(IN) :: n_phi_points !! Number of steps taken in angular direction in one of six triangles that make up the hexagon

  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :, :, :) !! Supercondicting pairings read from a simulation
  REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :) !! Charge density read from a simulation

  REAL(REAL64) :: Kappa_nearest(3) !! K-variables aligned with orbital's directions for nearest neighbours
  REAL(REAL64) :: Kappa_next(3) !! K-variables aligned with orbital's directions for next nearest neighbours
  COMPLEX(REAL64) :: C_l(3) !! Shape functions for each orbital
  COMPLEX(REAL64) :: Gamma_nearest_orb(3) !! Fourier transformed Gamma for given orbital at current k-point for nearest neighbour
  COMPLEX(REAL64) :: Gamma_next_orb(3) !! Fourier transformed Gamma for given orbital at current k-point for next nearest neighbour
  INTEGER(INT64), PARAMETER :: N_PROJECTIONS = 18 !! Number of distict projections in orbital \otimes spatial basis
  COMPLEX(REAL64) :: Projections_k_space_nearest(N_PROJECTIONS) = CMPLX(0.0d0, 0.0d0, KIND=REAL64) !! Accumulator for integrated projections over k-space for nearest neighbours
  COMPLEX(REAL64) :: Projections_k_space_next(N_PROJECTIONS) = CMPLX(0.0d0, 0.0d0, KIND=REAL64) !! Accumulator for integrated projections over k-space for next nearest neighbours
  COMPLEX(REAL64) :: Projections_real_space(N_PROJECTIONS) = CMPLX(0.0d0, 0.0d0, KIND=REAL64) !! Projections onto irreducible representations from real-space Gammas
  COMPLEX(REAL64) :: Gamma_flat(N_PROJECTIONS)!! Flattened array of gammas in real space, to be used in real-space resolved symmetries
  COMPLEX(REAL64) :: Active_orbital(3) !! Orbital for which to write projection. Used to utilize projection callbacks in writing basis functions for each orbital

  !File operations
  CHARACTER(LEN=200) :: filename
  INTEGER(INT32), PARAMETER:: FIRST_FILE_UNIT = 10
  INTEGER(INT32) :: file_count, gamma_nearest_weighted_file, gamma_next_weighted_file
  CHARACTER(LEN=4) :: Projections_name_mapping(N_PROJECTIONS)

  !K-space variables
  INTEGER(INT32) :: n_triangle, i_r, j_phi
  REAL(REAL64) :: phi_k, r_k, r_max, dr
  REAL(REAL64) :: kx, ky

  !Iterators
  INTEGER(INT32) :: n, neigh, orb, spin, layer, band

  ! Interface for basis functions
  ABSTRACT INTERFACE
    FUNCTION projection_interface(K_orb, Gamma_projected) RESULT(projection)
      USE, INTRINSIC :: iso_fortran_env, ONLY: real64
      REAL(REAL64), INTENT(IN) :: K_orb(3)
      COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(3)
      COMPLEX(REAL64) :: projection
    END FUNCTION projection_interface
  END INTERFACE
  ! Define a pointer to function type
  TYPE projection_function_cb_t
    PROCEDURE(projection_interface), POINTER, NOPASS :: cb
  END TYPE projection_function_cb_t
  ! Array of basis function pointers
  TYPE(projection_function_cb_t) :: Projections_cb(N_PROJECTIONS)

  CALL ASSIGN_PROJECTION_CALLBACKS(Projections_cb)
  CALL ASSIGN_PROJECTION_NAMES(Projections_name_mapping)

  CALL GET_INPUT(TRIM(projections % path)//"input.nml", sc_input)

  !Redefining step in radial and angular directions according to postprocessing parameters
  sc_input % discretization % derived % dr_k = R_K_MAX / projections % Nr_points
  sc_input % discretization % derived % dphi_k = (PI / 3.0d0) / projections % Nphi_points  !Slicing every hexagon's triangle into the same number of phi steps

  ALLOCATE (Gamma_SC(sc_input % discretization % ORBITALS, &
                    & N_ALL_NEIGHBOURS, &
                    & SPINS, &
                    & SPINS, &
                    & sc_input % discretization % derived % LAYER_COUPLINGS, &
                    & sc_input % discretization % SUBBANDS))
  ALLOCATE (Charge_dens(sc_input % discretization % derived % DIM_POSITIVE_K, &
                       & sc_input % discretization % SUBBANDS))

  CALL GET_SAFE_GAMMA_SC(Gamma_SC, projections % path, sc_input % discretization)
  CALL GET_SAFE_CHARGE_DENS(Charge_dens, projections % path, sc_input % discretization)

  ! Gamma_SC = 0.0d0
  ! Gamma_SC(:, :, 1, :, :) = 1 * meV2au
  ! Gamma_SC(:, :, 2, :, :) = -1 * meV2au

  CALL GET_REAL_SPACE_PROJECTIONS(Projections_real_space, Gamma_SC)
  DO n = 1, N_PROJECTIONS
    PRINT *, "Projection real space ", TRIM(ADJUSTL(Projections_name_mapping(n))), " = ", Projections_real_space(n) / SUM(Projections_real_space)
  END DO

  !Open files for projections
  file_count = FIRST_FILE_UNIT
  !Open files for projection basis functions
  DO n = 1, N_PROJECTIONS
    WRITE (filename, "(A, A)") "BasisFunc_", TRIM(ADJUSTL(Projections_name_mapping(n)))
    OPEN (unit=file_count, FILE=TRIM(projections % path)//"OutputData/"//TRIM(filename)//".dat", FORM="FORMATTED", ACTION="WRITE")
    WRITE (file_count, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_basis_orb_nearest_1)[meV]   Im(Gamma_basis_orb_nearest_1)[meV]   Re(Gamma_basis_orb_next_1)[meV]   Im(Gamma_basis_orb_next_1)[meV]  Re(Gamma_basis_orb_nearest_2)[meV] ... Re(Gamma_basis_weighted)[meV]   Im(Gamma_basis_weighted)[meV]"
    file_count = file_count + 1
  END DO
  gamma_nearest_weighted_file = file_count
  file_count = file_count + 1
  OPEN (unit=gamma_nearest_weighted_file, FILE=TRIM(projections % path)//"OutputData/Gamma_K_nearest_weighted.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_nearest_weighted_file, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_orb_1)[meV]   Im(Gamma_orb_1)[meV]  Re(Gamma_basis_2)[meV] ... Re(Gamma_weighted)[meV]   Im(Gamma_weighted)[meV]"

  gamma_next_weighted_file = file_count
  OPEN (unit=gamma_next_weighted_file, FILE=TRIM(projections % path)//"OutputData/Gamma_K_next_weighted.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_next_weighted_file, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_orb_1)[meV]   Im(Gamma_orb_1)[meV]  Re(Gamma_basis_2)[meV] ... Re(Gamma_weighted)[meV]   Im(Gamma_weighted)[meV]"

  !$omp parallel do collapse(3) schedule(dynamic, 1) private(phi_k, r_k, r_max, dr, kx, ky, Kappa_nearest, Kappa_next, C_l, Gamma_nearest_orb, Gamma_next_orb, n_triangle, j_phi, i_r, orb, n, Active_orbital)
  DO n_triangle = -N_BZ_SECTIONS / 2, N_BZ_SECTIONS / 2 - 1
    DO j_phi = 0, projections % Nphi_points - 1
      DO i_r = 0, projections % Nr_points
        !Must stuck everhing here to provide pure loops for collapse(3)
        phi_k = n_triangle * (PI / 3.0d0) + j_phi * sc_input % discretization % derived % dphi_k
        r_max = r_max_phi(MOD(ABS(phi_k), PI / 3))
        dr = r_max / projections % Nr_points
        r_k = i_r * dr

        !Transform from graphene reciprocal lattice to kx and ky
        kx = r_k * COS(phi_k)
        ky = r_k * SIN(phi_k)

        !Calculate k-space coordinates aligned with orbitals' directions
        Kappa_nearest(1) = SQRT(3.0d0) / 2.0d0 * kx + ky / 2.0d0
        Kappa_nearest(2) = -SQRT(3.0d0) / 2.0d0 * kx + ky / 2.0d0
        Kappa_nearest(3) = -ky

        Kappa_next(1) = -SQRT(3.0d0) / 2.0d0 * kx + 3.0d0 * ky / 2.0d0
        Kappa_next(2) = -SQRT(3.0d0) / 2.0d0 * kx - 3.0d0 * ky / 2.0d0
        Kappa_next(3) = SQRT(3.0d0) * kx

        !Calculate shape functions for orbitals
        C_l(1) = (1 + Kappa_nearest(2) * Kappa_nearest(3)) / (1 + kx**2 + ky**2)
        C_l(2) = (1 + Kappa_nearest(1) * Kappa_nearest(3)) / (1 + kx**2 + ky**2)
        C_l(3) = (1 + Kappa_nearest(1) * Kappa_nearest(2)) / (1 + kx**2 + ky**2)

        Gamma_nearest_orb = CMPLX(0.0d0, 0.0d0, KIND=REAL64)
        Gamma_next_orb = CMPLX(0.0d0, 0.0d0, KIND=REAL64)
        !Performing simple Fourier transform of Gamma in each orbital
        DO orb = 1, sc_input % discretization % ORBITALS
          !Nearest neighbours
          Gamma_nearest_orb(orb) = Gamma_SC(orb, 1, 1, 2, 1, 1) * pairing_1(ky) + Gamma_SC(orb, 1, 1, 2, 2, 1) * CONJG(pairing_1(ky)) + &
                                   Gamma_SC(orb, 2, 1, 2, 1, 1) * pairing_2(kx, ky) + Gamma_SC(orb, 2, 1, 2, 2, 1) * CONJG(pairing_2(kx, ky)) + &
                                   Gamma_SC(orb, 3, 1, 2, 1, 1) * pairing_3(kx, ky) + Gamma_SC(orb, 3, 1, 2, 2, 1) * CONJG(pairing_3(kx, ky))

          !Next-nearest neighbours
          Gamma_next_orb(orb) = Gamma_SC(orb, N_NEIGHBOURS + 1, 1, 2, 1, 1) * pairing_nnn_1(kx) + Gamma_SC(orb, N_NEIGHBOURS + 2, 1, 2, 1, 1) * pairing_nnn_2(kx, ky) + &
                                Gamma_SC(orb, N_NEIGHBOURS + 3, 1, 2, 1, 1) * pairing_nnn_3(kx, ky) + Gamma_SC(orb, N_NEIGHBOURS + 4, 1, 2, 1, 1) * pairing_nnn_4(kx) + &
                                Gamma_SC(orb, N_NEIGHBOURS + 5, 1, 2, 1, 1) * pairing_nnn_5(kx, ky) + Gamma_SC(orb, N_NEIGHBOURS + 6, 1, 2, 1, 1) * pairing_nnn_6(kx, ky)
        END DO

        !Write Gamma(k) resulting from the file that was read
        WRITE (gamma_nearest_weighted_file, '(10F15.5)') kx, ky, (REAL(Gamma_nearest_orb(orb)), AIMAG(Gamma_nearest_orb(orb)), orb=1, sc_input % discretization % ORBITALS),&
        & REAL(SUM(C_l * Gamma_nearest_orb)), AIMAG(SUM(C_l * Gamma_nearest_orb))

        WRITE (gamma_next_weighted_file, '(10F15.5)') kx, ky, (REAL(Gamma_next_orb(orb)), AIMAG(Gamma_next_orb(orb)), orb=1, sc_input % discretization % ORBITALS),&
        & REAL(SUM(C_l * Gamma_next_orb)), AIMAG(SUM(C_l * Gamma_next_orb))

        !Project onto basis functions and write basis functions onto which we project
        DO n = 1, N_PROJECTIONS
          !$omp critical (projections_update)
          Projections_k_space_nearest(n) = Projections_k_space_nearest(n) + Projections_cb(n) % cb(Kappa_nearest, Gamma_nearest_orb) * r_k * dr * sc_input % discretization % derived % dphi_k
          Projections_k_space_next(n) = Projections_k_space_next(n) + Projections_cb(n) % cb(Kappa_next, Gamma_next_orb) * r_k * dr * sc_input % discretization % derived % dphi_k
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') kx, ky
          DO orb = 1, sc_input % discretization % ORBITALS
            Active_orbital = 0.0d0
            Active_orbital(orb) = 1.0d0
            WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_nearest, Active_orbital)), AIMAG(Projections_cb(n) % cb(Kappa_nearest, Active_orbital))
            WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_next, Active_orbital)), AIMAG(Projections_cb(n) % cb(Kappa_next, Active_orbital))
          END DO
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_nearest, C_l)), AIMAG(Projections_cb(n) % cb(Kappa_nearest, C_l))
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)') REAL(Projections_cb(n) % cb(Kappa_next, C_l)), AIMAG(Projections_cb(n) % cb(Kappa_next, C_l))
        END DO
        !$omp end critical (projections_update)
      END DO
    END DO
  END DO
  !$omp end parallel do

  DO n = 1, N_PROJECTIONS
    PRINT *, "Projection k-space ", TRIM(ADJUSTL(Projections_name_mapping(n))), " = ", Projections_k_space_nearest(n) / SUM(Projections_k_space_nearest)
  END DO

  !Writing normalized projections since we are interested only in relative contributions
  OPEN (unit=gamma_next_weighted_file + 1, FILE=TRIM(projections % path)//"OutputData/Projections_nearest.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_next_weighted_file + 1, '(A)') "#Irrep      ABS(R-space projection)      ABS(K-space projection)"
  DO n = 1, N_PROJECTIONS
    WRITE (gamma_next_weighted_file + 1, '(A, 2F15.5)') Projections_name_mapping(n), &
    & ABS(Projections_real_space(n) / SUM(Projections_real_space)), ABS(Projections_k_space_nearest(n) / SUM(Projections_k_space_nearest)), &
    & ABS(Projections_k_space_next(n) / SUM(Projections_k_space_next))
  END DO
  CLOSE (gamma_next_weighted_file + 1)

  CLOSE (gamma_nearest_weighted_file)
  CLOSE (gamma_next_weighted_file)
  DO n = 1, N_PROJECTIONS
    CLOSE (FIRST_FILE_UNIT + n - 1)
  END DO

  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)

! Internal functions
CONTAINS

  SUBROUTINE ASSIGN_PROJECTION_CALLBACKS(Projections_callbacks)
    !! Assign callbacks to projection onto irreducible representations functions
    IMPLICIT NONE
    TYPE(projection_function_cb_t), INTENT(OUT) :: Projections_callbacks(N_PROJECTIONS) !! Array of projection callbacks

    !Assigning callbacks of projections
    Projections_callbacks(1) % cb => a1_1_proj
    Projections_callbacks(2) % cb => a1_2_proj
    Projections_callbacks(3) % cb => a2_1_proj
    Projections_callbacks(4) % cb => b1_1_proj
    Projections_callbacks(5) % cb => b1_2_proj
    Projections_callbacks(6) % cb => b2_1_proj
    Projections_callbacks(7) % cb => e1_1_proj
    Projections_callbacks(8) % cb => e1_2_proj
    Projections_callbacks(9) % cb => e1_3_proj
    Projections_callbacks(10) % cb => e1_4_proj
    Projections_callbacks(11) % cb => e1_5_proj
    Projections_callbacks(12) % cb => e1_6_proj
    Projections_callbacks(13) % cb => e2_1_proj
    Projections_callbacks(14) % cb => e2_2_proj
    Projections_callbacks(15) % cb => e2_3_proj
    Projections_callbacks(16) % cb => e2_4_proj
    Projections_callbacks(17) % cb => e2_5_proj
    Projections_callbacks(18) % cb => e2_6_proj
  END SUBROUTINE ASSIGN_PROJECTION_CALLBACKS

  SUBROUTINE ASSIGN_PROJECTION_NAMES(Projections_names)
    !! Map index of projection to its name
    IMPLICIT NONE
    CHARACTER(LEN=4) :: Projections_names(N_PROJECTIONS) !! Array of projection names mapped to projection number

    Projections_names(1) = "A1_1"
    Projections_names(2) = "A1_2"
    Projections_names(3) = "A2_1"
    Projections_names(4) = "B1_1"
    Projections_names(5) = "B1_2"
    Projections_names(6) = "B2_1"
    Projections_names(7) = "E1_1"
    Projections_names(8) = "E1_2"
    Projections_names(9) = "E1_3"
    Projections_names(10) = "E1_4"
    Projections_names(11) = "E1_5"
    Projections_names(12) = "E1_6"
    Projections_names(13) = "E2_1"
    Projections_names(14) = "E2_2"
    Projections_names(15) = "E2_3"
    Projections_names(16) = "E2_4"
    Projections_names(17) = "E2_5"
    Projections_names(18) = "E2_6"

  END SUBROUTINE ASSIGN_PROJECTION_NAMES

  SUBROUTINE GET_REAL_SPACE_PROJECTIONS(Projections_real_space, Gamma_SC)
    IMPLICIT NONE
    COMPLEX(REAL64), INTENT(OUT) :: Projections_real_space(N_PROJECTIONS) !! Array of projections in real space
    COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(sc_input % discretization % ORBITALS, &
                                      & N_ALL_NEIGHBOURS, &
                                      & SPINS, &
                                      & SPINS, &
                                      & sc_input % discretization % derived % LAYER_COUPLINGS, &
                                      & sc_input % discretization % SUBBANDS)
    COMPLEX(REAL64) :: Gamma_flat(N_PROJECTIONS)
    INTEGER(INT32) :: orb, neigh, n
    n = 1
    DO orb = 1, sc_input % discretization % ORBITALS
      DO neigh = 1, N_NEIGHBOURS
        Gamma_flat(n) = Gamma_SC(orb, neigh, 1, 2, 1, 1)
        n = n + 1
        !Get neighbor from opposite sublattice corresponding to the next vector in conter-clockwise rotation
        Gamma_flat(n) = Gamma_SC(orb, MOD(neigh + 1, N_NEIGHBOURS) + 1, 1, 2, 2, 1)
        n = n + 1
      END DO
    END DO
    Projections_real_space(1) = (Gamma_flat(3) + Gamma_flat(6) + Gamma_flat(8) + Gamma_flat(11) + Gamma_flat(13) + Gamma_flat(16))
    Projections_real_space(2) = (Gamma_flat(1) + Gamma_flat(2) + Gamma_flat(4) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(9) + Gamma_flat(10) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(15) + Gamma_flat(17) + Gamma_flat(18))
    Projections_real_space(3) = (Gamma_flat(2) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(10) + Gamma_flat(15) + Gamma_flat(18) - (Gamma_flat(1) + Gamma_flat(4) + Gamma_flat(9) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(17)))
    Projections_real_space(4) = (Gamma_flat(6) + Gamma_flat(8) + Gamma_flat(16) - (Gamma_flat(3) + Gamma_flat(11) + Gamma_flat(13)))
    Projections_real_space(5) = (Gamma_flat(2) + Gamma_flat(4) + Gamma_flat(10) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(18) - (Gamma_flat(1) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(9) + Gamma_flat(15) + Gamma_flat(17)))
    Projections_real_space(6) = (Gamma_flat(1) + Gamma_flat(2) + Gamma_flat(9) + Gamma_flat(10) + Gamma_flat(17) + Gamma_flat(18) - (Gamma_flat(4) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(15)))
    Projections_real_space(7) = (Gamma_flat(5) + Gamma_flat(10) - (Gamma_flat(2) + Gamma_flat(7)))
    Projections_real_space(8) = (Gamma_flat(6) + Gamma_flat(11) - (Gamma_flat(3) + Gamma_flat(8)))
    Projections_real_space(9) = (Gamma_flat(1) + Gamma_flat(12) - (Gamma_flat(4) + Gamma_flat(9)))
    Projections_real_space(10) = (Gamma_flat(3) + Gamma_flat(16) - (Gamma_flat(6) + Gamma_flat(13)))
    Projections_real_space(11) = (Gamma_flat(4) + Gamma_flat(17) - (Gamma_flat(1) + Gamma_flat(14)))
    Projections_real_space(12) = (Gamma_flat(5) + Gamma_flat(18) - (Gamma_flat(2) + Gamma_flat(15)))
    Projections_real_space(13) = (Gamma_flat(7) + Gamma_flat(10) - Gamma_flat(2) - Gamma_flat(5))
    Projections_real_space(14) = (Gamma_flat(8) + Gamma_flat(11) - Gamma_flat(3) - Gamma_flat(6))
    Projections_real_space(15) = (Gamma_flat(9) + Gamma_flat(12) - Gamma_flat(1) - Gamma_flat(4))
    Projections_real_space(16) = (Gamma_flat(13) + Gamma_flat(16) - Gamma_flat(3) - Gamma_flat(6))
    Projections_real_space(17) = (Gamma_flat(14) + Gamma_flat(17) - Gamma_flat(1) - Gamma_flat(4))
    Projections_real_space(18) = (Gamma_flat(15) + Gamma_flat(18) - Gamma_flat(2) - Gamma_flat(5))

    !Normalize to per-bond-per-orbital coupling in meV
    Projections_real_space = Projections_real_space / N_PROJECTIONS

  END SUBROUTINE GET_REAL_SPACE_PROJECTIONS

  !dir$ attributes forceinline :: a1_1_proj
  PURE FUNCTION a1_1_proj(k_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto A1_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(3) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * COS(K_orb(1))
    F_orb(2) = 2 * COS(K_orb(2))
    F_orb(3) = 2 * COS(K_orb(3))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION a1_1_proj

  !dir$ attributes forceinline :: a1_2_proj
  PURE FUNCTION a1_2_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto A1_2 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS)
    F_orb(1) = 2 * (COS(K_orb(2)) + COS(K_orb(3)))
    F_orb(2) = 2 * (COS(K_orb(1)) + COS(K_orb(3)))
    F_orb(3) = 4 * COS(K_orb(3) / 2.0d0) * COS((K_orb(1) - K_orb(2)) / 2.0d0)
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION a1_2_proj

  !dir$ attributes forceinline :: a2_1_proj
  PURE FUNCTION a2_1_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto A2_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * (COS(K_orb(2)) - COS(K_orb(3)))
    F_orb(2) = 2 * (-COS(K_orb(1)) + COS(K_orb(3)))
    F_orb(3) = 4 * SIN(K_orb(3) / 2.0d0) * SIN((K_orb(1) - K_orb(2)) / 2.0d0)
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION a2_1_proj

  !dir$ attributes forceinline :: b1_1_proj
  PURE FUNCTION b1_1_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto B1_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(1))
    F_orb(2) = 2 * imag * SIN(K_orb(2))
    F_orb(3) = 2 * imag * SIN(K_orb(3))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION b1_1_proj

  !dir$ attributes forceinline :: b1_2_proj
  PURE FUNCTION b1_2_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto B1_2 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * (SIN(K_orb(2)) + SIN(K_orb(3)))
    F_orb(2) = 2 * imag * (SIN(K_orb(1)) + SIN(K_orb(3)))
    F_orb(3) = -4 * imag * SIN(K_orb(3) / 2.0d0) * COS((K_orb(1) - K_orb(2)) / 2.0d0)
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION b1_2_proj

  !dir$ attributes forceinline :: b2_1_proj
  PURE FUNCTION b2_1_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto B2_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * (SIN(K_orb(2)) - SIN(K_orb(3)))
    F_orb(2) = 2 * imag * (-SIN(K_orb(1)) + SIN(K_orb(3)))
    F_orb(3) = 4 * imag * COS(K_orb(3) / 2.0d0) * SIN((K_orb(1) - K_orb(2)) / 2.0d0)
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION b2_1_proj

  !dir$ attributes forceinline :: e1_1_proj
  PURE FUNCTION e1_1_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(2))
    F_orb(2) = 2 * imag * SIN(K_orb(3))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_1_proj

  !dir$ attributes forceinline :: e1_2_proj
  PURE FUNCTION e1_2_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_2 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(1))
    F_orb(2) = -2 * imag * SIN(K_orb(2))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_2_proj

  !dir$ attributes forceinline :: e1_3_proj
  PURE FUNCTION e1_3_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_3 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(3))
    F_orb(2) = 2 * imag * SIN(K_orb(1))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_3_proj

  !dir$ attributes forceinline :: e1_4_proj
  PURE FUNCTION e1_4_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_4 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(1))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * imag * SIN(K_orb(3))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_4_proj

  !dir$ attributes forceinline :: e1_5_proj
  PURE FUNCTION e1_5_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_5 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(3))
    F_orb(2) = 0.d0
    F_orb(3) = -2 * imag * SIN(K_orb(2))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_5_proj

  !dir$ attributes forceinline :: e1_6_proj
  PURE FUNCTION e1_6_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E1_6 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(2))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * imag * SIN(K_orb(1))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e1_6_proj

  !dir$ attributes forceinline :: e2_1_proj
  PURE FUNCTION e2_1_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_1 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(2))
    F_orb(2) = 2 * COS(K_orb(3))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_1_proj

  !dir$ attributes forceinline :: e2_2_proj
  PURE FUNCTION e2_2_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_2 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(1))
    F_orb(2) = 2 * COS(K_orb(2))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_2_proj

  !dir$ attributes forceinline :: e2_3_proj
  PURE FUNCTION e2_3_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_3 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(3))
    F_orb(2) = 2 * COS(K_orb(1))
    F_orb(3) = 0.d0
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_3_proj

  !dir$ attributes forceinline :: e2_4_proj
  PURE FUNCTION e2_4_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_4 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(1))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(3))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_4_proj

  !dir$ attributes forceinline :: e2_5_proj
  PURE FUNCTION e2_5_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_5 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(3))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(2))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_5_proj

  !dir$ attributes forceinline :: e2_6_proj
  PURE FUNCTION e2_6_proj(K_orb, Gamma_projected) RESULT(proj)
    !! Project Gamma_projected onto E2_6 irreducible representation
    IMPLICIT NONE
    COMPLEX(REAL64) :: proj
    REAL(REAL64), INTENT(IN) :: K_orb(sc_input % discretization % ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX(REAL64), INTENT(IN) :: Gamma_projected(sc_input % discretization % ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX(REAL64) :: F_orb(sc_input % discretization % ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(2))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(1))
    proj = DOT_PRODUCT(F_orb, Gamma_projected)
  END FUNCTION e2_6_proj

END SUBROUTINE CALCULATE_PROJECTIONS
END MODULE symmetry
