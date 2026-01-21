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

PROGRAM MAIN
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE logger
USE hamiltonians
USE parameters
USE utilities
USE writers
USE reader
USE broydenV2
USE local_integrand
USE integrate
USE self_consistency
USE interaction_factory
USE types
USE omp_lib

IMPLICIT NONE

COMPLEX(REAL64), ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :), U_transformation(:, :)
REAL(REAL64), ALLOCATABLE :: Energies(:)
#ifndef BAND_BASIS
COMPLEX(REAL64), ALLOCATABLE :: Delta_new(:, :, :, :)
COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :, :), Gamma_SC_new(:, :, :, :)
#else
COMPLEX(REAL64), ALLOCATABLE :: Delta_new(:, :, :)
COMPLEX(REAL64), ALLOCATABLE :: Gamma_SC(:, :, :), Gamma_SC_new(:, :, :)
#endif
REAL(REAL64), ALLOCATABLE :: Delta_broyden(:), Delta_new_broyden(:)
REAL(REAL64), ALLOCATABLE :: Charge_dens(:, :), Charge_dens_new(:, :)

TYPE(sc_input_params_t) :: sc_input

REAL(REAL64) :: gamma_max_error, charge_max_error
REAL(REAL64) :: gamma_max_error_prev, charge_max_error_prev

INTEGER(INT32) :: band
INTEGER(INT32) :: sc_iter
LOGICAL :: sc_flag

INTEGER(INT32) :: delta_real_elems

!OMP specific
INTEGER(INT32) :: max_num_threads

max_num_threads = omp_get_max_threads()

CALL INIT_LOGGER()
WRITE (log_string, *) "Max num threads", max_num_threads
LOG_INFO(log_string)

CALL GET_INPUT("./input.nml", sc_input)
ASSOCIATE (SUBLATTICES => sc_input % discretization % SUBLATTICES, &
         & SUBBANDS => sc_input % discretization % SUBBANDS, &
         & ORBITALS => sc_input % discretization % ORBITALS, &
         & TBA_DIM => sc_input % discretization % derived % TBA_DIM, &
         & DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
         & DIM => sc_input % discretization % derived % DIM, &
         & LAYER_COUPLINGS => sc_input % discretization % derived % LAYER_COUPLINGS)
  !Basis
  !c_{k,yz,Ti1,up}, c_{k,zx,Ti1,up}, c_{k,xy,Ti1,up},
  !c_{k,yz,Ti2,up}, c_{k,zx,Ti2,up}, c_{k,xy,Ti2,up},
  !c_{k,yz,Ti1,down}, c_{k,zx,Ti1,down}, c_{k,xy,Ti1,down},
  !c_{k,yz,Ti2,down}, c_{k,zx,Ti2,down}, c_{k,xy,Ti2,down},
  ! + H.c{-k}
  ALLOCATE (Hamiltonian(DIM, DIM))
  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (Hamiltonian_const_band(DIM, DIM))
  ALLOCATE (U_transformation(DIM, DIM))
  ALLOCATE (Energies(DIM))
#ifndef BAND_BASIS
  ALLOCATE (Delta_new(N_ALL_NEIGHBOURS + N_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Gamma_SC(N_ALL_NEIGHBOURS + N_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Gamma_SC_new(N_ALL_NEIGHBOURS + N_NEIGHBOURS, DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#else
  ALLOCATE (Delta_new(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Gamma_SC(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Gamma_SC_new(DIM_POSITIVE_K, DIM_POSITIVE_K, SUBBANDS))
#endif
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Charge_dens_new(DIM_POSITIVE_K, SUBBANDS))
  !Fourth dimension for coupling between sublattices/layers
  !Coupling with nearest neighbours is inter-layer, thus we include both
  !Ti1 - Ti2 coupling and Ti2 - Ti1 coupling separately.
  !For next-to-nearest neighbours we only include Ti1-Ti1 etc. coupling
  !Due to its intra-layer character
  delta_real_elems = 2 * size(Gamma_SC) + size(Charge_dens)

  ALLOCATE (Delta_broyden(delta_real_elems))   !Flattened Gamma array
  ALLOCATE (Delta_new_broyden(delta_real_elems))
END ASSOCIATE

!Initializations
Hamiltonian = CMPLX(0., 0., KIND=REAL64)
Hamiltonian_const = CMPLX(0., 0., KIND=REAL64)
Hamiltonian_const_band = CMPLX(0., 0., KIND=REAL64)
U_transformation = CMPLX(0., 0., KIND=REAL64)
Energies = 0.

Delta_new = CMPLX(0., 0., KIND=REAL64)

ASSOCIATE (sc => sc_input % self_consistency, &
          & sb => sc_input % physical % subband_params)
  IF (sc % read_gamma_from_file) THEN
    LOG_INFO("Reading gamma from file: "//TRIM(sc % path_to_gamma_start))
    CALL GET_GAMMA_SC(Gamma_SC, TRIM(sc % path_to_gamma_start), sc_input % discretization)
  ELSE
    LOG_INFO("Initializing Gamma_SC")
    Gamma_SC = CMPLX(0., 0., KIND=REAL64)
    ! TODO: Think about different initialization techniques I can use here.
    CALL SET_GAMMA_INITIAL(Gamma_SC, &
                          & sb % J_tensor, &
                          & sc % gamma_start, &
                          & sc % gamma_nnn_start, &
                          & sc_input % discretization)
  END IF

  Gamma_SC_new = CMPLX(0., 0., KIND=REAL64)

  IF (sc % read_charge_from_file) THEN
    LOG_INFO("Reading charge from file: "//TRIM(sc % path_to_charge_start))
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(sc % path_to_charge_start), sc_input % discretization)
  ELSE
    Charge_dens = sc % charge_start
  END IF
END ASSOCIATE
Charge_dens_new = 0.

Delta_broyden = 0.
Delta_new_broyden = 0.

gamma_max_error = 0.
charge_max_error = 0.

!Computing k-independent terms
CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

OPEN (unit=99, FILE="./OutputData/Convergence.dat", FORM="FORMATTED", ACTION="WRITE")
WRITE (99, *) '# sc_iter, gamma_max_err, charge_max_err'
DO sc_iter = 1, sc_input % self_consistency % max_sc_iter
  WRITE (log_string, '(a, I0)') "==== SC_ITER: ", sc_iter
  LOG_INFO(log_string)
  DO band = 1, sc_input % discretization % SUBBANDS
    WRITE (log_string, '(a, I0)') " **** BAND: ", band
    LOG_INFO(log_string)

    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band, sc_input % physical % subband_params % Subband_energies, sc_input % discretization)

    !Integration over chunks is computed via Romberg algorithm.
#ifndef BAND_BASIS
    CALL ROMBERG_INTEGRATE_OVER_CHUNKS(Hamiltonian_const_band, Gamma_SC(:, :, :, band), Charge_dens(:, band), &
                                      & Delta_new(:, :, :, band), Charge_dens_new(:, band), sc_input)
#else
    CALL ROMBERG_INTEGRATE_OVER_CHUNKS(Hamiltonian_const_band, Gamma_SC(:, :, band), Charge_dens(:, band), &
                                      & Delta_new(:, :, band), Charge_dens_new(:, band), sc_input)
#endif
  END DO
  !Multiplying the result by Brillouin zone area
  Delta_new = Delta_new / JACOBIAN
  Charge_dens_new = Charge_dens_new / JACOBIAN
  !#########################################################################################################################
  !This is a critical section - only one thread can execute that and all thread should have ended their job up to that point
  !#########################################################################################################################
  !TODO: Check, because probably this is not needed anymor
  !CALL GET_GAMMAS_FROM_DELTAS(Gamma_SC_new, Delta_new, sc_input % discretization, sc_input % physical % subband_params % nearest_interorb_multiplier, sc_input % physical % subband_params % next_interorb_multiplier)
  Gamma_SC_new = Delta_new

  CALL COMPUTE_COOPER_PAIR_HOPPINGS(sc_input % physical % subband_params % cooper_pair_hoppings, &
    & sc_input % discretization % n_cooper_pair_hoppings, &
    & Gamma_SC_new, &
    & sc_input % discretization % derived % DIM_POSITIVE_K)

  CALL PRINT_GAMMA(Gamma_SC_new, "Gamma_SC_new", sc_input % discretization)

  CALL CHECK_CONVERGENCE(sc_flag, Gamma_SC, Gamma_SC_new, Charge_dens, Charge_dens_new, &
    & gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error, sc_iter, sc_input % self_consistency, sc_input % discretization)

  IF (sc_flag) THEN
    LOG_INFO("Convergence reached!")
    EXIT
  END IF

  WRITE (99, '(I0, 2E15.5)') sc_iter, gamma_max_error / meV2au, charge_max_error
  WRITE (log_string, '(a, E15.5)') "gamma_max_error [meV]: ", gamma_max_error / meV2au
  LOG_INFO(log_string)
  WRITE (log_string, '(a, E15.5)') "charge_max_error: ", charge_max_error
  LOG_INFO(log_string)

  !PRINT*, "Gamma max error ", gamma_max_error

  !In the beginning of convergence use Broyden method to quickly find minimum
  ! IF (gamma_max_error  > 1e-3) THEN
  !Broyden mixing
  !Flatten arrays for Broyden mixing
  CALL FLATTEN_FOR_BROYDEN(Gamma_SC, Gamma_SC_new, Charge_dens, Charge_dens_new, Delta_new_broyden, Delta_broyden, delta_real_elems, sc_input % discretization)

  CALL mix_broyden(delta_real_elems, Delta_new_broyden, Delta_broyden, sc_input % self_consistency % sc_alpha, sc_iter, 4, .FALSE.)

  CALL RESHAPE_FROM_BROYDEN(Gamma_SC, Charge_dens, Delta_broyden, delta_real_elems, sc_input % discretization)
  !In the last phase of convergence use linear mixing to avoid spare oscillations
  ! ELSE
  !     PRINT*, "Linear mixing"
  !Linear mixing
  ! Gamma_SC = (1.-sc_input % self_consistency % sc_alpha) * Gamma_SC + sc_input % self_consistency % sc_alpha * Gamma_SC_new
  ! Charge_dens = (1.-sc_input % self_consistency % sc_alpha) * Charge_dens + sc_input % self_consistency % sc_alpha * Charge_dens_new
  ! END IF

  !Gamma_SC(:,:,:,2) = CONJG(Gamma_SC(:,:,:,1)) !This is valid in asbence of magnetic field
  Delta_new = CMPLX(0., 0., KIND=REAL64)
  Gamma_SC_new = CMPLX(0., 0., KIND=REAL64)
  Charge_dens_new = 0.

  !To check the state of the simulation
  CALL PRINT_GAMMA(Gamma_SC, "Gamma_SC_iter", sc_input % discretization)
  CALL PRINT_CHARGE(Charge_dens, "Charge_dens_iter", sc_input % discretization)

END DO !End of SC loop
CLOSE (99)

CALL PRINT_GAMMA(Gamma_SC, "Gamma_SC_final", sc_input % discretization)
CALL PRINT_CHARGE(Charge_dens, "Charge_dens_final", sc_input % discretization)

!Just for memory deallocation, the .TRUE. flag is crucial
CALL mix_broyden(delta_real_elems, Delta_new_broyden, Delta_broyden, sc_input % self_consistency % sc_alpha, sc_iter, 4, .TRUE.)

CALL CLOSE_LOGGER()

DEALLOCATE (Hamiltonian)
DEALLOCATE (Hamiltonian_const)
DEALLOCATE (Hamiltonian_const_band)
DEALLOCATE (Energies)
DEALLOCATE (Delta_new)
DEALLOCATE (Gamma_SC)
DEALLOCATE (Gamma_SC_new)
DEALLOCATE (U_transformation)
DEALLOCATE (Delta_broyden)
DEALLOCATE (Delta_new_broyden)
DEALLOCATE (Charge_dens)
DEALLOCATE (Charge_dens_new)

ASSOCIATE (params => sc_input % physical % subband_params)
  IF (ALLOCATED(params % V_layer)) DEALLOCATE (params % V_layer) !Deallocate dynamic variable
  IF (ALLOCATED(params % Subband_energies)) DEALLOCATE (params % Subband_energies) !Deallocate dynamic variable
END ASSOCIATE

END PROGRAM MAIN
