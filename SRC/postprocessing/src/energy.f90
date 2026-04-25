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
MODULE energy
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

SUBROUTINE CALCULATE_DOS(dos_params)
  TYPE(post_dos_t), INTENT(IN) :: dos_params

  TYPE(sc_input_params_t) :: sc_input

  REAL(REAL64) :: E0
  INTEGER(INT32) :: DOS_steps
  INTEGER(INT32) :: hamiltonian_dim

  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :)
  REAL(REAL64), ALLOCATABLE :: Energies(:)

#ifndef BAND_BASIS
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :)
#else
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :)
#endif

  REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :)

  REAL(REAL64), ALLOCATABLE :: DOS(:, :), DOS_local(:, :)
  CHARACTER(LEN=20) :: output_format

  REAL(REAL64) :: k1, k2, kx, ky, dk1, dk2
  REAL(REAL64) :: sc_multiplier
  INTEGER(INT32) :: i, j, k, n, lat, orb, orb_prime, spin, band

  INTEGER(INT32) :: i_ref, j_ref
  REAL(REAL64) :: dk1_ref, dk2_ref

  INTEGER(INT32) :: points_within_energy_range, points_within_energy_range_local

  LOGICAL :: fileExists

  CALL GET_INPUT(TRIM(dos_params % path)//"input.nml", sc_input)

  DOS_steps = INT((dos_params % E_max - dos_params % E_min) / dos_params % dE0)
  dk1 = K1_MAX / dos_params % Nk_points
  dk2 = K2_MAX / dos_params % Nk_points

  dk1_ref = dk1 / (dos_params % Nk_points_refined + 1)
  dk2_ref = dk2 / (dos_params % Nk_points_refined + 1)

  ASSOCIATE (SUBLATTICES => sc_input % discretization % SUBLATTICES, &
         & SUBBANDS => sc_input % discretization % SUBBANDS, &
         & ORBITALS => sc_input % discretization % ORBITALS, &
         & TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
         & DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
         & DIM => sc_input % discretization % derived % DIM, &
         & LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)

    !If superconductivity is to be included, we add Nambu space to the Hamiltonian and double the size.
    IF (dos_params % include_sc) THEN
      hamiltonian_dim = DIM
      sc_multiplier = 0.5
    ELSE
      hamiltonian_dim = DIM_POSITIVE_K
      sc_multiplier = 1.0
    END IF

    ALLOCATE (Hamiltonian(DIM, DIM))
    ALLOCATE (Hamiltonian_const(DIM, DIM))
    ALLOCATE (Hamiltonian_const_band(DIM, DIM))
    ALLOCATE (Energies(hamiltonian_dim))
#ifndef BAND_BASIS
    ALLOCATE (Gamma_SC(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#else
    ALLOCATE (Gamma_SC(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#endif
    ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (DOS(hamiltonian_dim, 0:DOS_steps))
    ALLOCATE (DOS_local(hamiltonian_dim, 0:DOS_steps))
  END ASSOCIATE
  Hamiltonian = CMPLX(0., 0., KIND=REAL64)
  Hamiltonian_const = CMPLX(0., 0., KIND=REAL64)
  Energies = 0.
  Gamma_SC = CMPLX(0., 0., KIND=REAL64)
  Charge_dens = 0.
  DOS = 0.0d0
  DOS_local = 0.0d0

  CALL GET_SAFE_CHARGE_DENS(Charge_dens, dos_params % path, sc_input % discretization)

  IF (dos_params % include_sc) CALL GET_SAFE_GAMMA_SC(Gamma_SC, dos_params % path, sc_input % discretization)

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

  DO band = 1, sc_input % discretization % SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    !Adapt potential of given subband (energy difference due to quantization)
    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band, sc_input % physical % subband_params % Subband_energies, sc_input % discretization)

    WRITE (log_string, *) "Calculating energies and integrating DOS..."
    LOG_INFO(log_string)
    !$omp parallel private(E0, k1, k2, kx, ky, Hamiltonian, Energies, &
    !$omp                & DOS_local, points_within_energy_range_local)

    DOS_local = 0.0d0 !Initialize DOS_local for each thread!
    points_within_energy_range_local = 0
    !$omp do collapse(2)
    DO i = -dos_params % Nk_points / 2, dos_params % Nk_points / 2
      DO j = -dos_params % Nk_points / 2, dos_params % Nk_points / 2
        k1 = i * dk1
        k2 = j * dk2

        kx = 2.*PI / (SQRT(3.0d0)) * k1
        ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
        Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
        CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
        Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
        CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
        CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
        Hamiltonian = sc_multiplier * Hamiltonian !Should by multiplied by 0.5 if in Nambu space
        CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), hamiltonian_dim)

        ! If lowest energy is beyond the range we are calculating the DOS for, skip
        ! This might be changed if magnetic field is to be introduced
        IF (MINVAL(ABS(Energies)) > dos_params % E_max) CYCLE
        points_within_energy_range_local = points_within_energy_range_local + 1

        !Update DOS for current thread.
        DO n = 0, DOS_steps
          E0 = dos_params % E_min + n * dos_params % dE0
          DO k = 1, hamiltonian_dim
            DOS_local(k, n) = DOS_local(k, n) + dirac_delta(Energies(k), E0, dos_params % zeta_DOS)
          END DO
        END DO

        !Add a grid refinement here, if lowest energy is in [E_DOS_min, E_DOS_max]
        !Center of the cell
        IF (i .lt. dos_params % Nk_points / 2 .AND. j .lt. dos_params % Nk_points / 2) THEN
          DO i_ref = 1, dos_params % Nk_points_refined
            DO j_ref = 1, dos_params % Nk_points_refined
              k1 = i * dk1 + i_ref * dk1_ref
              k2 = j * dk2 + j_ref * dk2_ref
              kx = 2.*PI / (SQRT(3.0d0)) * k1
              ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
              Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
              CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
              Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
              CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
              CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
              Hamiltonian = sc_multiplier * Hamiltonian !Should by multiplied by 0.5 if in Nambu space
              CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), hamiltonian_dim)

              !Update DOS for current thread.
              DO n = 0, DOS_steps
                E0 = dos_params % E_min + n * dos_params % dE0
                DO k = 1, hamiltonian_dim
                  DOS_local(k, n) = DOS_local(k, n) + dirac_delta(Energies(k), E0, dos_params % zeta_DOS)
                END DO
              END DO

            END DO
          END DO
        END IF

        IF ((i .lt. dos_params % Nk_points / 2 .AND. j .lt. dos_params % Nk_points / 2) .OR. (i .eq. dos_params % Nk_points / 2)) THEN
          !Left edge without corner
          DO j_ref = 1, dos_params % Nk_points_refined
            k1 = i * dk1
            k2 = j * dk2 + j_ref * dk2_ref
            kx = 2.*PI / (SQRT(3.0d0)) * k1
            ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
            Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
            CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
            Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
            CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
            CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
            Hamiltonian = sc_multiplier * Hamiltonian !Should by multiplied by 0.5 if in Nambu space
            CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), hamiltonian_dim)

            !Update DOS for current thread.
            DO n = 0, DOS_steps
              E0 = dos_params % E_min + n * dos_params % dE0
              DO k = 1, hamiltonian_dim
                DOS_local(k, n) = DOS_local(k, n) + dirac_delta(Energies(k), E0, dos_params % zeta_DOS)
              END DO
            END DO

          END DO
        END IF

        IF ((i .lt. dos_params % Nk_points / 2 .AND. j .lt. dos_params % Nk_points / 2) .OR. (j .eq. dos_params % Nk_points / 2)) THEN
          !Bottom edge without corner
          DO i_ref = 1, dos_params % Nk_points_refined
            k1 = i * dk1 + i_ref * dk1_ref
            k2 = j * dk2
            kx = 2.*PI / (SQRT(3.0d0)) * k1
            ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
            Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
            CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
            Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
            CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
            CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
            Hamiltonian = sc_multiplier * Hamiltonian !Should by multiplied by 0.5 if in Nambu space
            CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), hamiltonian_dim)

            !Update DOS for current thread.
            DO n = 0, DOS_steps
              E0 = dos_params % E_min + n * dos_params % dE0
              DO k = 1, hamiltonian_dim
                DOS_local(k, n) = DOS_local(k, n) + dirac_delta(Energies(k), E0, dos_params % zeta_DOS)
              END DO
            END DO

          END DO
        END IF

      END DO
    END DO
    !$omp end do

    !$omp critical (accumulate_dos)
    DO n = 0, DOS_steps
      DOS(:, n) = DOS(:, n) + DOS_local(:, n)
    END DO
    points_within_energy_range = points_within_energy_range + points_within_energy_range_local
    !$omp end critical (accumulate_dos)

    !$omp end parallel
  END DO

  ! Normalize the DOS
  IF (MAXVAL(SUM(DOS, DIM=1)) .NE. 0.0d0) DOS = DOS / MAXVAL(SUM(DOS, DIM=1))

  WRITE (log_string, *) "Included points: ", points_within_energy_range, "out of ", dos_params % Nk_points**2
  LOG_INFO(log_string)

  WRITE (log_string, *) "Writing DOS to file"
  LOG_INFO(log_string)

  output_format = '(2E15.5, *(E15.5))'  !Adjust according to number of bands
  OPEN (unit=9, FILE=TRIM(dos_params % path)//"OutputData/DOS.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, '(A)') "#E[meV]   DOS_total[a.u]    DOS_band1[a.u]   DOS_band2[a.u]   ... DOS_bandN[a.u]"
  DO n = 0, DOS_steps
    E0 = dos_params % E_min + n * dos_params % dE0
    WRITE (9, output_format) E0 / meV2au, SUM(DOS(:, n)), (DOS(k, n), k=1, hamiltonian_dim)
  END DO
  CLOSE (9)

  DEALLOCATE (Hamiltonian)
  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (Hamiltonian_const_band)
  DEALLOCATE (Energies)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (DOS)
  DEALLOCATE (DOS_local)
  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_DOS

SUBROUTINE CALCULATE_DISPERSION(dispersion)
  !! Calculates dispersion relation in the first Brillouin zone.
  !! Takes physical parameters from input.nml from a directory specified by inputPath.
  TYPE(post_dispersion_relation_t), INTENT(INOUT) :: dispersion
  ! CHARACTER(LEN=*), INTENT(IN) :: inputPath
  ! INTEGER(INT32), INTENT(IN) :: Nk_points
  ! LOGICAL, INTENT(IN) :: include_sc

  TYPE(sc_input_params_t) :: sc_input
  CHARACTER(LEN=20) :: output_format

  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :)
  COMPLEX(REAL64), ALLOCATABLE :: U_k(:, :), U_minus_k(:, :), U_combined(:, :), Hamiltonian_dummy(:, :)
  REAL(REAL64), ALLOCATABLE :: Energies(:)

#ifndef BAND_BASIS
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :)
#else
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :)
#endif
  REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :)

  REAL(REAL64) :: kx, ky
  REAL(REAL64) :: phi_k, r_max, dr, r_k
  REAL(REAL64) :: sc_multiplier
  INTEGER(INT32) :: i, j, k, n, lat, orb, orb_prime, spin, l, m, band
  INTEGER(INT32) :: i_r, j_phi, n_triangle
  INTEGER(INT32) :: s_up_idx, s_down_idx
  INTEGER(INT32) :: kx_steps, ky_steps
  REAL(REAL64) :: yz_contribution, zx_contribution, xy_contribution
  REAL(REAL64), ALLOCATABLE :: Lat_contributions(:)
  REAL(REAL64) :: spin_x_contribution, spin_y_contribution, spin_z_contribution
  REAL(REAL64) :: electron_contribution, hole_contribution

  LOGICAL :: fileExists
  INTEGER(INT32) :: hamiltonian_dim

  CALL GET_INPUT(TRIM(dispersion % path)//"input.nml", sc_input)

  !Redefining step in radial and angular directions according to postprocessing parameters
  sc_input % discretization % derived % dr_k = R_K_MAX / dispersion % Nr_points
  sc_input % discretization % derived % dphi_k = (PI / 3.0d0) / dispersion % Nphi_points  !Slicing every hexagon's triangle into the same number of phi steps

  !If superconductivity is to be included, we add Nambu space to the Hamiltonian and double the size.
  IF (dispersion % include_sc) THEN
    hamiltonian_dim = sc_input % discretization % derived % DIM
    sc_multiplier = 0.5d0
  ELSE
    hamiltonian_dim = sc_input % discretization % derived % DIM_POSITIVE_K
    sc_multiplier = 1.0d0
  END IF

  ASSOCIATE (SUBLATTICES => sc_input % discretization % SUBLATTICES, &
        & SUBBANDS => sc_input % discretization % SUBBANDS, &
        & ORBITALS => sc_input % discretization % ORBITALS, &
        & TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
        & DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
        & DIM => sc_input % discretization % derived % DIM, &
        & LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)

    WRITE (output_format, '(A, I0, A)') '(I5, ', 11 + SUBLATTICES, 'E15.5)'

    ALLOCATE (Hamiltonian(DIM, DIM))
    ALLOCATE (Hamiltonian_const(DIM, DIM))
    ALLOCATE (Hamiltonian_const_band(DIM, DIM))
    ALLOCATE (Hamiltonian_dummy(DIM, DIM))
    ALLOCATE (Energies(hamiltonian_dim))
#ifndef BAND_BASIS
    ALLOCATE (Gamma_SC(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#else
    ALLOCATE (Gamma_SC(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#endif
    ALLOCATE (U_k(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (U_minus_k(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (U_combined(DIM, DIM))
    ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
    ALLOCATE (Lat_contributions(SUBLATTICES))
  END ASSOCIATE

  Hamiltonian = CMPLX(0., 0., KIND=REAL64)
  Hamiltonian_const = CMPLX(0., 0., KIND=REAL64)
  Energies = 0.
  Gamma_SC = CMPLX(0., 0., KIND=REAL64) * meV2au
  Charge_dens = 0.
  Lat_contributions = 0.

  IF ((sc_input % physical % subband_params % U_HUB .NE. 0.0) .OR. &
     & (sc_input % physical % subband_params % V_HUB .NE. 0.0)) THEN
    CALL GET_SAFE_CHARGE_DENS(Charge_dens, dispersion % path, sc_input % discretization)
  END IF

  IF (dispersion % include_sc) THEN
    CALL GET_SAFE_GAMMA_SC(Gamma_SC, dispersion % path, sc_input % discretization)
  END IF

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

  OPEN (unit=9, FILE=TRIM(dispersion % path)//"OutputData/Energies.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, '(A)') "#N kx[1/a] ky[1/a] Energy[meV] P(yz) P(zx) P(xy) P(lat1) P(lat2) ... P(latN) P(sigma_x) P(sigma_y) P(sigma_z) P(electron) P(hole)"

  DO band = 1, sc_input % discretization % SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    !Adapt potential of given subband (energy difference due to quantization)
    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band, sc_input % physical % subband_params % Subband_energies, sc_input % discretization)

    !$omp parallel do collapse(3) schedule(dynamic, 1) private(phi_k, r_k, r_max, dr, kx, ky, n_triangle, j_phi, i_r, orb, n, &
    !$omp                                                    & Energies, Hamiltonian, yz_contribution, zx_contribution, xy_contribution, &
    !$omp                                                    & Lat_contributions, spin_x_contribution, spin_y_contribution, spin_z_contribution, &
    !$omp                                                    & electron_contribution, hole_contribution, s_up_idx, s_down_idx, l, m, spin, lat, &
    !$omp                                                    & U_k, U_minus_k, U_combined, Hamiltonian_dummy)
    DO n_triangle = -N_BZ_SECTIONS / 2, N_BZ_SECTIONS / 2 - 1
      DO j_phi = 0, dispersion % Nphi_points - 1
        DO i_r = 0, dispersion % Nr_points
          phi_k = n_triangle * (PI / 3.0d0) + j_phi * sc_input % discretization % derived % dphi_k
          r_max = r_max_phi(MOD(ABS(phi_k), PI / 3))
          dr = r_max / dispersion % Nr_points
          r_k = i_r * dr

          !Transform from graphene reciprocal lattice to kx and ky
          kx = r_k * COS(phi_k)
          ky = r_k * SIN(phi_k)

          Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
          Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
          CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
          U_k = Hamiltonian(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                          & :sc_input % discretization % derived % DIM_POSITIVE_K)
          U_minus_k = Hamiltonian(sc_input % discretization % derived % DIM_POSITIVE_K + 1:, &
                                & sc_input % discretization % derived % DIM_POSITIVE_K + 1:)

          CALL DIAGONALIZE_HERMITIAN(U_k, &
                                   & Energies(:sc_input % discretization % derived % DIM_POSITIVE_K), &
                                   & sc_input % discretization % derived % DIM_POSITIVE_K)
          CALL DIAGONALIZE_HERMITIAN(U_minus_k, &
                                   & Energies(:sc_input % discretization % derived % DIM_POSITIVE_K), &
                                   & sc_input % discretization % derived % DIM_POSITIVE_K)
          U_combined = CMPLX(0., 0., KIND=REAL64)
          U_combined(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                   & :sc_input % discretization % derived % DIM_POSITIVE_K) = U_k
          U_combined(sc_input % discretization % derived % DIM_POSITIVE_K + 1:, &
                   & sc_input % discretization % derived % DIM_POSITIVE_K + 1:) = CONJG(U_minus_k)

          CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
          Hamiltonian = sc_multiplier * Hamiltonian !Should by multiplied by 0.5 if in Nambu space

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), hamiltonian_dim)

#ifdef BAND_BASIS
          Hamiltonian_dummy = Hamiltonian
          Hamiltonian(:hamiltonian_dim, :hamiltonian_dim) = MATMUL(U_combined(:hamiltonian_dim, :hamiltonian_dim), &
            & Hamiltonian_dummy(:hamiltonian_dim, :hamiltonian_dim))
#endif

          !Calculate contributions
          DO l = 1, hamiltonian_dim
            !Distinguishing orbital contributions
            yz_contribution = 0.
            zx_contribution = 0.
            xy_contribution = 0.
            DO n = 1, hamiltonian_dim, sc_input % discretization % ORBITALS
              yz_contribution = yz_contribution + ABS(Hamiltonian(n, l))**2
              zx_contribution = zx_contribution + ABS(Hamiltonian(n + 1, l))**2
              xy_contribution = xy_contribution + ABS(Hamiltonian(n + 2, l))**2
            END DO

            !Distinguishing lattice contributions
            Lat_contributions(:) = 0.
            DO m = 0, sc_input % discretization % SUBLATTICES - 1
              DO spin = 0, 1
                DO n = 1, sc_input % discretization % ORBITALS
                  Lat_contributions(m + 1) = Lat_contributions(m + 1) + &
                  & ABS(Hamiltonian(spin * sc_input % discretization % derived % TBA_DIM + m * sc_input % discretization % ORBITALS + n, l))**2
                  IF (dispersion % include_sc) THEN
                    Lat_contributions(m + 1) = Lat_contributions(m + 1) + &
                    & ABS(Hamiltonian(sc_input % discretization % derived % DIM_POSITIVE_K + spin * sc_input % discretization % derived % TBA_DIM + m * sc_input % discretization % ORBITALS + n, l))**2
                  END IF
                END DO
              END DO
            END DO

            !Distinguishing spin contributions
            spin_x_contribution = 0.
            spin_y_contribution = 0.
            spin_z_contribution = 0.
            DO n = 1, sc_input % discretization % derived % TBA_DIM
              s_up_idx = n
              s_down_idx = n + sc_input % discretization % derived % TBA_DIM
              !sigma_x
              spin_x_contribution = spin_x_contribution + CONJG(Hamiltonian(s_up_idx, l)) * Hamiltonian(s_down_idx, l)
              spin_x_contribution = spin_x_contribution + CONJG(Hamiltonian(s_down_idx, l)) * Hamiltonian(s_up_idx, l)
              !sigma_y
              spin_y_contribution = spin_y_contribution - imag * CONJG(Hamiltonian(s_up_idx, l)) * Hamiltonian(s_down_idx, l)
              spin_y_contribution = spin_y_contribution + imag * CONJG(Hamiltonian(s_down_idx, l)) * Hamiltonian(s_up_idx, l)
              !sigma_z
              spin_z_contribution = spin_z_contribution + ABS(Hamiltonian(s_up_idx, l))**2 !Spin up
              spin_z_contribution = spin_z_contribution - ABS(Hamiltonian(s_down_idx, l))**2 !Spin down
              IF (dispersion % include_sc) THEN
                s_up_idx = n + sc_input % discretization % derived % DIM_POSITIVE_K
                s_down_idx = n + sc_input % discretization % derived % TBA_DIM + sc_input % discretization % derived % DIM_POSITIVE_K
                !sigma_x
                spin_x_contribution = spin_x_contribution + CONJG(Hamiltonian(s_up_idx, l)) * Hamiltonian(s_down_idx, l)
                spin_x_contribution = spin_x_contribution + CONJG(Hamiltonian(s_down_idx, l)) * Hamiltonian(s_up_idx, l)
                !sigma_y
                spin_y_contribution = spin_y_contribution - imag * CONJG(Hamiltonian(s_up_idx, l)) * Hamiltonian(s_down_idx, l)
                spin_y_contribution = spin_y_contribution + imag * CONJG(Hamiltonian(s_down_idx, l)) * Hamiltonian(s_up_idx, l)
                !sigma_z
                spin_z_contribution = spin_z_contribution + ABS(Hamiltonian(s_up_idx, l))**2 !Spin up
                spin_z_contribution = spin_z_contribution - ABS(Hamiltonian(s_down_idx, l))**2 !Spin down
              END IF
            END DO

            electron_contribution = 0.
            hole_contribution = 0.
            IF (dispersion % include_sc) THEN
              DO n = 1, sc_input % discretization % derived % DIM_POSITIVE_K
                electron_contribution = electron_contribution + ABS(Hamiltonian(n, l))**2
                hole_contribution = hole_contribution + ABS(Hamiltonian(sc_input % discretization % derived % DIM_POSITIVE_K + n, l))**2
              END DO
            ELSE
              electron_contribution = 1.
              hole_contribution = 0.
            END IF

            WRITE (9, output_format) (band - 1) * hamiltonian_dim + l, kx, ky, Energies(l) / meV2au, &
            & yz_contribution, zx_contribution, xy_contribution, &
            & (Lat_contributions(lat), lat=1, sc_input % discretization % SUBLATTICES), &
            & spin_x_contribution, spin_y_contribution, spin_z_contribution, &
            & electron_contribution, hole_contribution
          END DO
        END DO
      END DO
    END DO
    !$omp end parallel do
  END DO
  CLOSE (9)

  DEALLOCATE (Hamiltonian)
  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (Energies)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (Lat_contributions)
  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_DISPERSION

SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP(gap)

  TYPE(post_sc_gap_t) :: gap
  ! CHARACTER(LEN=*), INTENT(IN) :: inputPath !! This should be a path to folder where input.nml resides
  ! REAL(REAL64), INTENT(IN) :: dE
  ! INTEGER(INT32), INTENT(IN) :: nBrillouinPoints
  TYPE(sc_input_params_t) :: sc_input
  CHARACTER(LEN=20) :: output_format

  COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :)
  REAL(REAL64), ALLOCATABLE :: Energies(:)

#ifndef BAND_BASIS
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :)
#else
  COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :)
#endif
  REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :)

  INTEGER(INT8), ALLOCATABLE :: IsFermiSurface(:, :) !! This indicates whether given (kx,ky) point is at Fermi surface
  INTEGER(INT16), ALLOCATABLE :: OrbitalAtFermiSurface(:, :) !! This indicates which state consitutes to the Fermi surface.
                                                         !! Because ZHEEV is used, energies are sorted from lowest to highest.
                                                         !! This implies that we do not recognize which spin, orbital etc.
                                                         !! constitutes to the Fermi surface, only "lowest", "second lowest" and so on.

  REAL(REAL64) :: brillouinZoneVertices(6, 2)

  REAL(REAL64) :: k1, k2, kx, ky, dkx, dky
  INTEGER(INT32) :: i, j, n, m, band
  INTEGER(INT32) :: kx_steps, ky_steps

  LOGICAL :: fileExists

  brillouinZoneVertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  brillouinZoneVertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  dkx = KX_MAX / gap % Nk_points
  dky = KY_MAX / gap % Nk_points

  kx_steps = INT(gap % Nk_points)
  ky_steps = INT(gap % Nk_points)

  !Get parameters from simulation
  CALL GET_INPUT(TRIM(gap % path)//"input.nml", sc_input)
  ASSOCIATE (SUBLATTICES => sc_input % discretization % SUBLATTICES, &
         & SUBBANDS => sc_input % discretization % SUBBANDS, &
         & ORBITALS => sc_input % discretization % ORBITALS, &
         & TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
         & DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
         & DIM => sc_input % discretization % derived % DIM, &
         & LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)

    ALLOCATE (Hamiltonian(DIM, DIM))
    ALLOCATE (Hamiltonian_const(DIM, DIM))
    ALLOCATE (Hamiltonian_const_band(DIM, DIM))
    ALLOCATE (Energies(DIM))
    ALLOCATE (IsFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
    ALLOCATE (OrbitalAtFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
#ifndef BAND_BASIS
    ALLOCATE (Gamma_SC(N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#else
    ALLOCATE (Gamma_SC(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#endif
    ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
  END ASSOCIATE
  Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
  Hamiltonian_const(:, :) = CMPLX(0., 0., KIND=REAL64)
  Energies(:) = 0.
  Gamma_SC = CMPLX(0., 0., KIND=REAL64) * meV2au
  Charge_dens = 0.

  output_format = '(3E15.5, I10)'

  !Calculation of superconducting gap
  CALL GET_SAFE_GAMMA_SC(Gamma_SC, gap % path, sc_input % discretization)

  CALL GET_SAFE_CHARGE_DENS(Charge_dens, gap % path, sc_input % discretization)
  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

  OPEN (unit=9, FILE=TRIM(gap % path)//"OutputData/SuperconductingGap.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) '#kx[1/a] ky[1/a] gap_SC[meV] N_orbital'
  DO band = 1, sc_input % discretization % SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    IsFermiSurface(:, :) = 0
    OrbitalAtFermiSurface(:, :) = 0

    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band, sc_input % physical % subband_params % Subband_energies, sc_input % discretization)

    !Dispersion relation in a normal state
    !$omp parallel private(kx, ky, Hamiltonian, Energies)
    !$omp do
    DO i = -kx_steps, kx_steps
      DO j = -ky_steps, ky_steps
        kx = i * dkx !* (2. * PI * 2./3.)
        ky = j * dky !* (2. * PI * 2./3.)

        IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

          Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
          Energies(:) = 0.
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
          Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
          CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
          CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:sc_input % discretization % derived % DIM_POSITIVE_K, &
                                    & :sc_input % discretization % derived % DIM_POSITIVE_K), &
                                    & Energies(:sc_input % discretization % derived % DIM_POSITIVE_K), sc_input % discretization % derived % DIM_POSITIVE_K)
          !Check whether current wavevector is in the Fermi surface
          IF (MINVAL(ABS(Energies(:sc_input % discretization % derived % DIM_POSITIVE_K))) < gap % dE) THEN
            IsFermiSurface(i, j) = 1
            OrbitalAtFermiSurface(i, j) = (band - 1) * sc_input % discretization % derived % DIM_POSITIVE_K + MINLOC(ABS(Energies(:sc_input % discretization % derived % DIM_POSITIVE_K)), 1)
          END IF
        END IF
      END DO
    END DO
    !$omp end do
    !$omp end parallel

    WRITE (log_string, *) "Fermi surface done"
    LOG_INFO(log_string)

    !$omp parallel private(kx, ky, Hamiltonian, Energies)
    !$omp do
    DO i = -kx_steps, kx_steps
      DO j = -ky_steps, ky_steps
        IF (IsFermiSurface(i, j) == 1) THEN
          kx = i * dkx
          ky = j * dky

          Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
          Energies(:) = 0.
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
          Hamiltonian = Hamiltonian_const_band + Hamiltonian
#ifndef BAND_BASIS
          CALL COMPUTE_INTERACTIONS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#else
          CALL COMPUTE_INTERACTIONS_BAND_BASIS(Hamiltonian, kx, ky, Charge_dens, Gamma_SC, sc_input % discretization, sc_input % physical)
#endif
          Hamiltonian = 0.5 * Hamiltonian !Should by multiplied by 0.5 if in Nambu space
          CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, sc_input % discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle
          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:, :), Energies(:), sc_input % discretization % derived % DIM)
          !Write superconducting gap
          WRITE (9, output_format) kx, ky, ABS(Energies(sc_input % discretization % derived % DIM_POSITIVE_K) - &
          & Energies(sc_input % discretization % derived % DIM_POSITIVE_K + 1)) / meV2au, OrbitalAtFermiSurface(i, j)
        END IF
      END DO
    END DO
    !$omp end do
    !$omp end parallel
  END DO
  CLOSE (9)

  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP

!########################### HELPER FUNCTIONS ###########################

LOGICAL FUNCTION is_inside_polygon(verticesArray, nVertices, pointX, pointY)
    !! This function returns true if the point (pointX, pointY) is inside the polygon
    !! Where vertices are defined in array verticesArray, containing X and Y coordinates.
    !! The assumption is that the polygon is two-dimensional
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: nVertices !! Number of vertices of polygon
  REAL(REAL64), INTENT(IN) :: verticesArray(nVertices, 2) !! X and Y coordinates of vertices
  REAL(REAL64), INTENT(IN) :: pointX, pointY !! X and Y coordinates of point to be tested
  INTEGER(INT32) :: i, j
  REAL(REAL64) :: xIntersection

  is_inside_polygon = .FALSE.

  DO i = 1, nVertices

    IF (pointX == verticesArray(i, 1) .AND. pointY == verticesArray(i, 2)) THEN
      is_inside_polygon = .TRUE.
      RETURN
    END IF

    j = MOD(i, nVertices) + 1 !! Index of next vertex
        !! Check if horizontal line of y = pointY can cross the line joining two vertices
    IF (pointY > MIN(verticesArray(i, 2), verticesArray(j, 2)) .AND. &
    & pointY <= MAX(verticesArray(i, 2), verticesArray(j, 2))) THEN
      IF (pointX <= MAX(verticesArray(i, 1), verticesArray(j, 1))) THEN
        IF (verticesArray(i, 2) /= verticesArray(j, 2)) THEN !To avoid division by zero
          !Calculates intersection point of line connecting two vertices and horizontal line y = pointY
          xIntersection = (pointY - verticesArray(i, 2)) * &
          & (verticesArray(j, 1) - verticesArray(i, 1)) / (verticesArray(j, 2) &
          & - verticesArray(i, 2)) + verticesArray(i, 1)

          IF (pointX <= xIntersection .OR. verticesArray(i, 1) == verticesArray(j, 1)) THEN
            is_inside_polygon = .NOT. is_inside_polygon
          END IF
        END IF
      END IF
    END IF
  END DO

  RETURN
END FUNCTION is_inside_polygon

END MODULE energy
