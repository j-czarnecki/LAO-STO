#include "macros_def.f90"

PROGRAM MAIN

USE mod_logger
USE mod_hamiltonians
USE mod_parameters
USE mod_utilities
USE mod_writers
USE mod_reader
USE mod_broydenV2
USE mod_compute_hamiltonians
USE mod_integrate
USE mod_self_consistency
USE omp_lib

IMPLICIT NONE

COMPLEX*16, ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :), U_transformation(:, :)
REAL*8, ALLOCATABLE :: Energies(:)
COMPLEX*16, ALLOCATABLE :: Delta_local(:, :, :, :, :), Delta_new(:, :, :, :, :)
REAL*8, ALLOCATABLE :: Delta_broyden(:), Delta_new_broyden(:)
COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :), Gamma_SC_new(:, :, :, :, :)
REAL*8, ALLOCATABLE :: Charge_dens(:, :), Charge_dens_new(:, :), Charge_dens_local(:, :)

REAL*8 :: gamma_error, gamma_max_error, charge_error, charge_max_error
REAL*8 :: gamma_max_error_prev, charge_max_error_prev

INTEGER*4 :: i, j, n, lat, orb, orb_prime, spin, band, band_prime
INTEGER*4 :: sc_iter
LOGICAL :: sc_flag

INTEGER*4 :: counter
INTEGER*4 :: delta_real_elems
INTEGER*4 :: broyden_index

!OMP specific
INTEGER*4 :: max_num_threads, used_threads

max_num_threads = omp_get_max_threads()

CALL INIT_LOGGER("")
WRITE (log_string, *) "Max num threads", max_num_threads
LOG_INFO(log_string)

! CALL omp_set_num_threads(max_num_threads)
! !$omp parallel
! PRINT*, "Hello from process", omp_get_thread_num()
! used_threads = omp_get_num_threads()
! CALL SLEEP(3)
! !$omp end parallel
! PRINT*, "Allocated threads", used_threads

CALL GET_INPUT("./input.nml")
!N_NEIGHBOURS + N_NEXT_NEIGHBOURS = 9, to implement both pairings
!2 because spin up-down and down-up,
!2 because of complex number
!LAYER_COUPLINGS because of Ti1-Ti2 coupling and Ti2 - Ti1 coupling (stored in this order) for nearest-neighbours
!DIM_POSITIVE_K included due to Charge density self-consistency
delta_real_elems = SUBBANDS * (DIM_POSITIVE_K + ORBITALS * 2 * 2 * (N_NEIGHBOURS * LAYER_COUPLINGS + N_NEXT_NEIGHBOURS * SUBLATTICES))

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
ALLOCATE (Delta_local(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))    !Third dimension for spin coupling up-down and down-up
!Fourth dimension for coupling between sublattices/layers
!Coupling with nearest neighbours is inter-layer, thus we include both
!Ti1 - Ti2 coupling and Ti2 - Ti1 coupling separately.
!For next-to-nearest neighbours we only include Ti1-Ti1 etc. coupling
!Due to its intra-layer character
ALLOCATE (Delta_new(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
ALLOCATE (Gamma_SC_new(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
ALLOCATE (Delta_broyden(delta_real_elems))   !Flattened Gamma array
ALLOCATE (Delta_new_broyden(delta_real_elems))
ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
ALLOCATE (Charge_dens_local(DIM_POSITIVE_K, SUBBANDS))
ALLOCATE (Charge_dens_new(DIM_POSITIVE_K, SUBBANDS))

!Initializations
Hamiltonian = DCMPLX(0., 0.)
Hamiltonian_const = DCMPLX(0., 0.)
Hamiltonian_const_band = DCMPLX(0., 0.)
U_transformation = DCMPLX(0., 0.)
Energies = 0.

Delta_local = DCMPLX(0., 0.)
Delta_new = DCMPLX(0., 0.)

IF (read_gamma_from_file) THEN
  LOG_INFO("Reading gamma from file: "//TRIM(path_to_gamma_start))
  CALL GET_GAMMA_SC(Gamma_SC, TRIM(path_to_gamma_start))
ELSE
  !Breaking spin up-down down-up symmetry
  !coupling for nearest-neighbours
  Gamma_SC(:, :N_NEIGHBOURS, 1, :, :) = DCMPLX(gamma_start, 0.)
  Gamma_SC(:, :N_NEIGHBOURS, 2, :, :) = DCMPLX(-gamma_start, 0.)
  !coupling for next nearest neighbours
  Gamma_SC(:, (N_NEIGHBOURS + 1):, 1, :, :) = DCMPLX(gamma_nnn_start, 0.)
  Gamma_SC(:, (N_NEIGHBOURS + 1):, 2, :, :) = DCMPLX(-gamma_nnn_start, 0.)
END IF

Gamma_SC_new = DCMPLX(0., 0.)

IF (read_charge_from_file) THEN
  LOG_INFO("Reading charge from file: "//TRIM(path_to_charge_start))
  CALL GET_CHARGE_DENS(Charge_dens, TRIM(path_to_charge_start))
ELSE
  Charge_dens = charge_start
END IF

Charge_dens_new = 0.
Charge_dens_local = 0.

Delta_broyden = 0.
Delta_new_broyden = 0.

gamma_max_error = 0.
charge_max_error = 0.

!Computing k-independent terms
CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:, :))
CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:, :))
CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:, :))
CALL COMPUTE_LAYER_POTENTIAL(Hamiltonian_const(:, :))
CALL COMPUTE_FERMI_ENERGY(Hamiltonian_const(:, :))
CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

OPEN (unit=99, FILE="./OutputData/Convergence.dat", FORM="FORMATTED", ACTION="WRITE")
WRITE (99, *) '# sc_iter, Re(Gamma), Im(Gamma), Re(Gamma_new), Im(Gamma_new), n, n_new, gamma_max_err, charge_max_err'
DO sc_iter = 1, max_sc_iter
  WRITE (log_string, '(a, I0)') "==== SC_ITER: ", sc_iter
  LOG_INFO(log_string)
  DO band = 1, SUBBANDS
    WRITE (log_string, '(a, I0)') " **** BAND: ", band
    LOG_INFO(log_string)

    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band)

    !Integration over chunks is computed via Romberg algorithm.
    !$omp parallel do collapse(2) schedule(dynamic, 1) private(Delta_local, Charge_dens_local)
    DO i = -k1_steps / 2, k1_steps / 2 - 1
      DO j = -k2_steps / 2, k2_steps / 2 - 1
        ! WRITE(log_string, *) 'Integrating over chunk: ', i, j
        ! LOG_INFO(log_string)

        CALL ROMBERG_Y(Hamiltonian_const_band(:, :), Gamma_SC(:, :, :, :, band), Charge_dens(:, band), i * dk1, (i + 1) * dk1, j * dk2, (j + 1) * dk2, &
        & Delta_local(:, :, :, :, band), Charge_dens_local(:, band), romb_eps_x, interpolation_deg_x, max_grid_refinements_x, &
        & romb_eps_y, interpolation_deg_y, max_grid_refinements_y)
        !This has to be atomic operations, since Delta_new and Charge_dens would be global variables for all threads
        !$omp critical (update_delta_and_charge)
        Delta_new(:, :, :, :, band) = Delta_new(:, :, :, :, band) + Delta_local(:, :, :, :, band)
        Charge_dens_new(:, band) = Charge_dens_new(:, band) + Charge_dens_local(:, band)
        !$omp end critical (update_delta_and_charge)
      END DO
    END DO !End of k-loop
    !$omp end parallel do
  END DO
  !#########################################################################################################################
  !This is a critical section - only one thread can execute that and all thread should have ended their job up to that point
  !#########################################################################################################################

  CALL GET_GAMMAS_FROM_DELTAS(Gamma_SC_new, Delta_new)
  CALL CHECK_CONVERGENCE(sc_flag, Gamma_SC, Gamma_SC_new, Charge_dens, Charge_dens_new, gamma_max_error_prev, gamma_max_error, charge_max_error_prev, charge_max_error, sc_iter)

  IF (sc_flag) THEN
    LOG_INFO("Convergence reached!")
    EXIT
  END IF

  WRITE (99, '(I0, 8E15.5)') sc_iter, REAL(Gamma_SC(1, 1, 1, 1, 1) / meV2au), AIMAG(Gamma_SC(1, 1, 1, 1, 1) / meV2au), &
  &                                 REAL(Gamma_SC_new(1, 1, 1, 1, 1) / meV2au), AIMAG(Gamma_SC_new(1, 1, 1, 1, 1) / meV2au), &
  &                                 Charge_dens(1, 1), Charge_dens_new(1, 1), gamma_max_error / meV2au, charge_max_error
  WRITE (log_string, '(a, E15.5)') "gamma_max_error [meV]: ", gamma_max_error / meV2au
  LOG_INFO(log_string)
  WRITE (log_string, '(a, E15.5)') "charge_max_error: ", charge_max_error
  LOG_INFO(log_string)

  !PRINT*, "Gamma max error ", gamma_max_error

  !In the beginning of convergence use Broyden method to quickly find minimum
  ! IF (gamma_max_error  > 1e-3) THEN
  !Broyden mixing
  !Flatten arrays for Broyden mixing
  CALL FLATTEN_FOR_BROYDEN(Gamma_SC, Gamma_SC_new, Charge_dens, Charge_dens_new, Delta_new_broyden, Delta_broyden, delta_real_elems)

  CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .FALSE.)

  CALL RESHAPE_FROM_BROYDEN(Gamma_SC, Charge_dens, Delta_broyden, delta_real_elems)
  !In the last phase of convergence use linear mixing to avoid spare oscillations
  ! ELSE
  !     PRINT*, "Linear mixing"
  !     !Linear mixing
  !     Gamma_SC(:,:,:,:) = (1. - sc_alpha)*Gamma_SC(:,:,:,:) + sc_alpha*Gamma_SC_new(:,:,:,:)
  ! END IF

  !Gamma_SC(:,:,:,2) = CONJG(Gamma_SC(:,:,:,1)) !This is valid in asbence of magnetic field
  Delta_new = DCMPLX(0., 0.)
  Gamma_SC_new = DCMPLX(0., 0.)
  Charge_dens_new = 0.

  !To check the state of the simulation
  CALL PRINT_GAMMA(Gamma_SC(:, :, :, :, :), "Gamma_SC_iter")
  CALL PRINT_CHARGE(Charge_dens(:, :), "Charge_dens_iter")

END DO !End of SC loop
CLOSE (99)

!Printing results after the simulation is done
CALL PRINT_GAMMA(Gamma_SC(:, :, :, :, :), "Gamma_SC_final")
CALL PRINT_CHARGE(Charge_dens(:, :), "Charge_dens_final")

!Just for memory deallocation, the .TRUE. flag is crucial
CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .TRUE.)

CALL CLOSE_LOGGER()

DEALLOCATE (Hamiltonian)
DEALLOCATE (Hamiltonian_const)
DEALLOCATE (Hamiltonian_const_band)
DEALLOCATE (Energies)
DEALLOCATE (Delta_local)
DEALLOCATE (Delta_new)
DEALLOCATE (Gamma_SC)
DEALLOCATE (Gamma_SC_new)
DEALLOCATE (U_transformation)
DEALLOCATE (Delta_broyden)
DEALLOCATE (Delta_new_broyden)
DEALLOCATE (Charge_dens)
DEALLOCATE (Charge_dens_local)
DEALLOCATE (Charge_dens_new)

IF (ALLOCATED(V_layer)) DEALLOCATE (V_layer) !Deallocate global variable
IF (ALLOCATED(Subband_energies)) DEALLOCATE (Subband_energies) !Deallocate global variable

END PROGRAM MAIN
