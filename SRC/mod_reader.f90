#include "macros_def.f90"

MODULE mod_reader
USE mod_parameters
USE mod_logger
IMPLICIT NONE
SAVE

!Default values, overwritten in get_input
!Discretization
INTEGER*4 :: k1_steps = 0
INTEGER*4 :: k2_steps = 0
INTEGER*4 :: SUBLATTICES = 2
INTEGER*4 :: SUBBANDS = 1

!Those parameters are in fact derived and recalculated if needed at GET_INPUT
!(see SET_HAMILTONIAN_PARAMS)
INTEGER*4 :: ORBITALS = 0
INTEGER*4 :: TBA_DIM = 0
INTEGER*4 :: DIM_POSITIVE_K = 0  !Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
INTEGER*4 :: DIM = 0    !*2 to transform to Nambu Space.
INTEGER*4 :: LAYER_COUPLINGS = 0 !This determines how many layer-related superconducting parameters have to be calculated.

!Physical parameters
REAL*8 :: T = 0.
REAL*8 :: t_D = 0.
REAL*8 :: t_I = 0.
REAL*8 :: t_Rashba = 0.
REAL*8 :: lambda_SOC = 0.
REAL*8 :: delta_trigonal = 0.
REAL*8 :: zeta_tetragonal = 0.
INTEGER*4 :: orb_affected_tetragonal = 1
REAL*8 :: v = 0.
REAL*8 :: V_pdp = 0.
REAL*8 :: V_pds = 0.
REAL*8 :: J_SC = 0.
REAL*8 :: J_SC_PRIME = 0.
REAL*8 :: J_SC_NNN = 0.
REAL*8 :: J_SC_PRIME_NNN = 0.
REAL*8 :: U_HUB = 0.
REAL*8 :: V_HUB = 0.
REAL*8 :: E_Fermi = 0.
REAL*8, ALLOCATABLE :: V_layer(:)
REAL*8, ALLOCATABLE :: Subband_energies(:)
REAL*8 :: B_field(3)

!Self-consistency
LOGICAL :: read_gamma_from_file = .FALSE.
CHARACTER(1000) :: path_to_gamma_start
LOGICAL :: read_charge_from_file = .FALSE.
CHARACTER(1000) :: path_to_charge_start
REAL*8 :: gamma_start = 0.
REAL*8 :: gamma_nnn_start = 0.
REAL*8 :: charge_start = 0.
INTEGER*4 :: max_sc_iter = 0
REAL*8 :: sc_alpha = 0.
REAL*8 :: sc_alpha_adapt = 0.
REAL*8 :: gamma_eps_convergence = 0.
REAL*8 :: charge_eps_convergence = 0.

!Romberg integration
REAL*8 :: romb_eps_x = 0.
INTEGER*4 :: interpolation_deg_x = 0
INTEGER*4 :: max_grid_refinements_x = 0
REAL*8 :: romb_eps_y = 0.
INTEGER*4 :: interpolation_deg_y = 0
INTEGER*4 :: max_grid_refinements_y = 0

!Derived
REAL * 8 eta_p
REAL*8 :: dr_k, dphi_k, domega

!Used for postprocessing
!Superconducting gap calculation
LOGICAL :: enable_sc_gap_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_sc_gap = ""
REAL*8 :: dE_sc_gap = 0.
INTEGER*4 :: Nk_points_sc_gap = 0
INTEGER*4 :: Nk_points_sc_gap_refined = 0

!For chern number calculation
LOGICAL :: enable_chern_number_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_chern_number = ""
INTEGER*4 :: Nk_points_chern_number = 0

!Dispersion relation calculation
LOGICAL :: enable_dispersion_relation_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_dispersion_relation = ""
LOGICAL :: include_sc_in_dispersion = .FALSE.
INTEGER*4 :: Nk_points_dispersion_relation = 0

!DOS calculation
LOGICAL :: enable_dos_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_dos = ""
REAL*8 :: E_DOS_min = 0
REAL*8 :: E_DOS_max = 0
REAL*8 :: dE0 = 0
REAL*8 :: zeta_DOS = 0
LOGICAL :: include_sc_in_dos = .FALSE.
INTEGER*4 :: Nk_points_dos = 0
INTEGER*4 :: Nk_points_dos_refined = 0

!Gamma as a map of k-vector calculation
LOGICAL :: enable_gamma_k_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_gamma_k = ""
INTEGER*4 :: Nk_points_gamma_k = 0

!Projections calculation
LOGICAL :: enable_projections_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_projections = ""
INTEGER*4 :: Nr_points_projections
INTEGER*4 :: Nphi_points_projections

NAMELIST /physical_params/  &
& T,                        &
& t_D,                      &
& t_I,                      &
& t_Rashba,                 &
& lambda_SOC,               &
& delta_trigonal,           &
& zeta_tetragonal,          &
& orb_affected_tetragonal,  &
& v,                        &
& V_pdp,                    &
& V_pds,                    &
& J_SC,                     &
& J_SC_PRIME,               &
& J_SC_NNN,                 &
& J_SC_PRIME_NNN,           &
& U_HUB,                    &
& V_HUB,                    &
& E_Fermi,                  &
& V_layer,                  &
& Subband_energies,         &
& B_field

NAMELIST /discretization/ &
& k1_steps,               &
& k2_steps,               &
& SUBLATTICES,            &
& SUBBANDS

NAMELIST /self_consistency/ &
& read_gamma_from_file,     &
& path_to_gamma_start,      &
& read_charge_from_file,    &
& path_to_charge_start,     &
& gamma_start,              &
& gamma_nnn_start,          &
& charge_start,             &
& max_sc_iter,              &
& sc_alpha,                 &
& sc_alpha_adapt,           &
& gamma_eps_convergence,    &
& charge_eps_convergence

NAMELIST /romberg_integration/ &
& romb_eps_x,                 &
& interpolation_deg_x,        &
& max_grid_refinements_x,     &
& romb_eps_y,                 &
& interpolation_deg_y,        &
& max_grid_refinements_y

NAMELIST /sc_gap_calculation/ &
& enable_sc_gap_calc,         &
& path_to_run_dir_sc_gap,     &
& dE_sc_gap,                  &
& Nk_points_sc_gap,           &
& Nk_points_sc_gap_refined

NAMELIST /chern_number_calculation/ &
& enable_chern_number_calc, &
& path_to_run_dir_chern_number, &
& Nk_points_chern_number

NAMELIST /dispersion_relation_calculation/ &
& enable_dispersion_relation_calc, &
& path_to_run_dir_dispersion_relation, &
& include_sc_in_dispersion, &
& Nk_points_dispersion_relation

NAMELIST /dos_calculation/ &
& enable_dos_calc, &
& path_to_run_dir_dos, &
& E_DOS_min, &
& E_DOS_max, &
& dE0, &
& zeta_DOS, &
& include_sc_in_dos, &
& Nk_points_dos, &
& Nk_points_dos_refined

NAMELIST /gamma_k_calculation/ &
& enable_gamma_k_calc, &
& path_to_run_dir_gamma_k, &
& Nk_points_gamma_k

NAMELIST /projections_calculation/ &
& enable_projections_calc, &
& path_to_run_dir_projections, &
& Nr_points_projections, &
& Nphi_points_projections

CONTAINS
SUBROUTINE GET_INPUT(nmlfile)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: nmlfile
  INTEGER*4 :: i
  INTEGER*4 :: io_status

  OPEN (unit=9, FILE=nmlfile, FORM="FORMATTED", ACTION="READ", STATUS="OLD")

  READ (9, NML=discretization, IOSTAT=io_status)
  IF (io_status .NE. 0) THEN
    WRITE (log_string, *) "Error reading discretization ", io_status
    LOG_ERROR(log_string)
    STOP "Error reading discretization"
  END IF
  WRITE (log_string, '(2(A, I0))') "k1_steps: ", k1_steps, &
                                & " k2_steps: ", k2_steps
  LOG_INFO(log_string)

  IF ((k1_steps .LE. 0) .OR. (k2_steps .LE. 0)) STOP "k_steps must be > 0"
  IF (SUBLATTICES .LE. 0) STOP "SUBLATTICES must be > 0"
  CALL SET_HAMILTONIAN_PARAMS()
  WRITE (log_string, '(6(A, I0))') " SUBLATTICES: ", SUBLATTICES,&
                                & " SUBBANDS: ", SUBBANDS,&
                                & " ORBITALS: ", ORBITALS, &
                                & " TBA_DIM: ", TBA_DIM, &
                                & " DIM_POSITIVE_K: ", DIM_POSITIVE_K, &
                                & " DIM: ", DIM
  LOG_INFO(log_string)

  !This is crucial
  ALLOCATE (V_layer(SUBLATTICES))
  ALLOCATE (Subband_energies(SUBBANDS))
  V_layer = 0.0d0
  Subband_energies = 0.0d0

  !TODO: WRITE BETTER CHECKS!!!!!!!!!!!!!!!!!!
  REWIND (9)
  READ (9, NML=physical_params, IOSTAT=io_status)
  IF (io_status .NE. 0) THEN
    WRITE (log_string, *) "Error reading physical_params ", io_status
    LOG_ERROR(log_string)
    STOP "Error reading physical_params"
  END IF

  WRITE (log_string, '(18(A, E15.5))') "T: ", T,&
                                   & " t_D: ", t_D,&
                                   & " t_I: ", t_I,&
                                   & " t_Rashba: ", t_Rashba,&
                                   & " lambda_SOC: ", lambda_SOC,&
                                   & " delta_trigonal: ", delta_trigonal,&
                                   & " zeta_tetragonal: ", zeta_tetragonal,&
                                   & " orb_affected_tetragonal: ", REAL(orb_affected_tetragonal),&
                                   & " v: ", v,&
                                   & " V_pdp: ", V_pdp,&
                                   & " V_pds: ", V_pds,&
                                   & " J_SC: ", J_SC,&
                                   & " J_SC_PRIME: ", J_SC_PRIME,&
                                   & " J_SC_NNN: ", J_SC_NNN,&
                                   & " J_SC_PRIME_NNN: ", J_SC_PRIME_NNN,&
                                   & " U_HUB: ", U_HUB,&
                                   & " V_HUB: ", V_HUB,&
                                   & " E_Fermi: ", E_Fermi
  LOG_INFO(log_string)
  WRITE (log_string, *) "V_layer: ", (V_layer(i), i=1, SUBLATTICES)
  LOG_INFO(log_string)
  WRITE (log_string, *) "Subband_energies: ", (Subband_energies(i), i=1, SUBBANDS)
  LOG_INFO(log_string)
  WRITE (log_string, '(A, 3E15.5)') "B_field: ", (B_field(i), i=1, 3)
  LOG_INFO(log_string)

  !Check input data
  IF (T < 0) STOP "Temperature in kelvins must be >= 0!"

  !Change to atomic units
  t_D = t_D * meV2au
  t_I = t_I * meV2au
  t_Rashba = t_Rashba * meV2au
  lambda_SOC = lambda_SOC * meV2au
  delta_trigonal = delta_trigonal * meV2au
  zeta_tetragonal = zeta_tetragonal * meV2au
  v = v * meV2au
  V_pdp = V_pdp * meV2au
  V_pds = V_pds * meV2au
  J_SC = J_SC * meV2au
  J_SC_PRIME = J_SC_PRIME * meV2au
  J_SC_NNN = J_SC_NNN * meV2au
  J_SC_PRIME_NNN = J_SC_PRIME_NNN * meV2au
  U_HUB = U_HUB * meV2au
  V_HUB = V_HUB * meV2au
  E_Fermi = E_Fermi * meV2au
  V_layer = V_layer * meV2au
  Subband_energies = Subband_energies * meV2au
  B_field = B_field * T2au

  REWIND (9)
  READ (9, NML=self_consistency, IOSTAT=io_status)
  IF (io_status .NE. 0) THEN
    WRITE (log_string, *) "Error reading self_consistency ", io_status
    LOG_ERROR(log_string)
    STOP "Error reading self_consistency"
  END IF
  IF (read_gamma_from_file .eqv. .FALSE.) WRITE (path_to_gamma_start, *) ""
  IF (read_charge_from_file .eqv. .FALSE.) WRITE (path_to_charge_start, *) ""

  WRITE (log_string, '(2(A, L1, 2A), 10(A, E15.5))') "read_gamma_from_file: ", read_gamma_from_file,&
                                                  & " path_to_gamma_start: ", TRIM(path_to_gamma_start),&
                                                  & " read_charge_from_file: ", read_charge_from_file,&
                                                  & " path_to_charge_start: ", TRIM(path_to_charge_start),&
                                                  & " gamma_start: ", gamma_start,&
                                                  & " gamma_nnn_start: ", gamma_nnn_start,&
                                                  & " charge_start: ", charge_start,&
                                                  & " max_sc_iter: ", REAL(max_sc_iter),&
                                                  & " sc_alpha: ", sc_alpha,&
                                                  & " sc_alpha_adapt: ", sc_alpha_adapt,&
                                                  & " gamma_eps_convergence: ", gamma_eps_convergence,&
                                                  & " charge_eps_convergence: ", charge_eps_convergence
  LOG_INFO(log_string)

  IF (read_gamma_from_file .AND. path_to_gamma_start == "") STOP "If read_gamma_from_file == .TRUE. then path_to_gamma_start must not be empty"
  IF (read_charge_from_file .AND. path_to_charge_start == "") STOP "If read_charge_from_file == .TRUE. then path_to_charge_start must not be empty"
  IF (charge_start .LT. 0) STOP "charge_start must be >= 0"
  IF (max_sc_iter .LE. 0) STOP "max_sc_iter must be > 0"
  IF (sc_alpha .LE. 0) STOP "sc_alpha (mixing parameter) must be > 0"
  IF ((sc_alpha_adapt .LE. 0) .OR. (sc_alpha_adapt .GT. 1)) STOP "sc_alpha_adapt must be in interval [0,1]"
  IF (gamma_eps_convergence .LE. 0) STOP "gamma_eps_convergence must be > 0"
  IF (charge_eps_convergence .LE. 0) STOP "charge_eps_convergence must be > 0"

  gamma_start = gamma_start * meV2au
  gamma_nnn_start = gamma_nnn_start * meV2au
  gamma_eps_convergence = gamma_eps_convergence * meV2au
  !Calculating derived values
  dr_k = R_K_MAX / k1_steps
  dphi_k = (PI / 3.0d0) / k2_steps  !Slicing every hexagon's triangle into the same number of phi steps
  !domega = ABS(dk1 * dk2 * SIN(2 * PI / 3.)) / (SIN(2 * PI / 3.) * K1_MAX * K2_MAX)
  eta_p = v * SQRT(3.) / 3.905 * nm2au

  REWIND (9)
  READ (9, NML=romberg_integration, IOSTAT=io_status)
  IF (io_status .NE. 0) THEN
    WRITE (log_string, *) "Error reading romberg_integration ", io_status
    LOG_ERROR(log_string)
    STOP "Error reading romberg_integration"
  END IF
  WRITE (log_string, '(2((A, E15.5), 2(A, I0)))') "romb_eps_x: ", romb_eps_x, &
                                               & " interpolation_deg_x: ", interpolation_deg_x, &
                                               & " max_grid_refinements_x: ", max_grid_refinements_x, &
                                               & " romb_eps_y: ", romb_eps_y, &
                                               & " interpolation_deg_y: ", interpolation_deg_y, &
                                               & " max_grid_refinements_y: ", max_grid_refinements_y
  LOG_INFO(log_string)

  IF (romb_eps_x .LE. 0) STOP "romb_eps_x must be > 0"
  IF (interpolation_deg_x .LE. 0) STOP "interpolation_deg_x must be > 0"
  IF (max_grid_refinements_x .LE. 0) STOP "max_grid_refinements_x must be > 0"

  IF (romb_eps_y .LE. 0) STOP "romb_eps_y must be > 0"
  IF (interpolation_deg_y .LE. 0) STOP "interpolation_deg_y must be > 0"
  IF (max_grid_refinements_y .LE. 0) STOP "max_grid_refinements_y must be > 0"

  CLOSE (9)

END SUBROUTINE GET_INPUT

SUBROUTINE GET_POSTPROCESSING_INPUT(nmlfile)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: nmlfile

  OPEN (unit=9, FILE=nmlfile, FORM="FORMATTED", ACTION="READ", STATUS="OLD")
  READ (9, NML=sc_gap_calculation)
  IF (enable_sc_gap_calc) THEN
    IF (Nk_points_sc_gap .LE. 0) STOP "Nk_points_sc_gap must be > 0"
    IF (path_to_run_dir_sc_gap == "") STOP "path_to_run_dir_sc_gap must not be empty"
    IF (dE_sc_gap .LE. 0) STOP "dE_sc_gap must be > 0"

    dE_sc_gap = dE_sc_gap * meV2au
  END IF

  READ (9, NML=chern_number_calculation)
  IF (enable_chern_number_calc) THEN
    IF (Nk_points_chern_number .LE. 0) STOP "Nk_points_chern_number must be > 0"
    IF (path_to_run_dir_chern_number == "") STOP "path_to_run_dir_chern_number must not be empty"
  END IF

  READ (9, NML=dispersion_relation_calculation)
  IF (enable_dispersion_relation_calc) THEN
    IF (Nk_points_dispersion_relation .LE. 0) STOP "Nk_points_dispersion_relation must be > 0"
    IF (path_to_run_dir_dispersion_relation == "") STOP "path_to_run_dir_dispersion_relation must not be empty"
  END IF

  READ (9, NML=dos_calculation)
  IF (enable_dos_calc) THEN
    IF (path_to_run_dir_dos == "") STOP "path_to_run_dir_dos must not be empty"
    IF (E_DOS_min .GT. E_DOS_max) STOP "E_DOS_min must be smaller than E_DOS_max"
    IF (dE0 .LE. 0) STOP "dE0 must be > 0"
    IF (zeta_DOS .LE. 0) STOP "zeta_DOS must be > 0"
    IF (Nk_points_dos .LE. 0) STOP "Nk_points_dos must be > 0"
    IF (Nk_points_dos_refined .LT. 0) STOP "Nk_points_dos_refined must be >= 0"
    E_DOS_min = E_DOS_min * meV2au
    E_DOS_max = E_DOS_max * meV2au
    zeta_DOS = zeta_DOS * meV2au
    dE0 = dE0 * meV2au
  END IF

  READ (9, NML=gamma_k_calculation)
  IF (enable_gamma_k_calc) THEN
    IF (Nk_points_gamma_k .LE. 0) STOP "Nk_points_gamma_k must be > 0"
    IF (path_to_run_dir_gamma_k == "") STOP "path_to_run_dir_gamma_k must not be empty"
  END IF

  READ (9, NML=projections_calculation)
  IF (enable_projections_calc) THEN
    IF (Nr_points_projections .LE. 0) STOP "Nr_points_projections must be > 0"
    IF (Nphi_points_projections .LE. 0) STOP "Nphi_points_projections must be > 0"
    IF (path_to_run_dir_projections == "") STOP "path_to_run_dir_projections must not be empty"
  END IF

  CLOSE (9)

END SUBROUTINE GET_POSTPROCESSING_INPUT

SUBROUTINE GET_SAFE_GAMMA_SC(Gamma_SC, input_path)
  !! This subroutine tries to read the Gamma_SC file from self-consistent simulation, checking for final result first,
  !! then checking the latest iteration file, and logging an error if neither is found.
  COMPLEX*16, INTENT(OUT) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS) !! Gamma array to be filled
  CHARACTER(LEN=*), INTENT(IN) :: input_path !! input path in which to look for the Gamma file
  LOGICAL :: file_exists !! whether the file exists or not

  !Read Gammas
  !Check the convergent simulation first
  INQUIRE (FILE=TRIM(input_path)//"OutputData/Gamma_SC_final.dat", EXIST=file_exists)
  IF (file_exists) THEN
    CALL GET_GAMMA_SC(Gamma_SC, TRIM(input_path)//"OutputData/Gamma_SC_final.dat")
  ELSE
    !Check the non-convergent simulation latest iteration
    INQUIRE (FILE=TRIM(input_path)//"OutputData/Gamma_SC_iter.dat", EXIST=file_exists)
    IF (file_exists) THEN
      CALL GET_GAMMA_SC(Gamma_SC, TRIM(input_path)//"OutputData/Gamma_SC_iter.dat")
    ELSE
      !No other file could have been produced, no Gamma file found
      WRITE (log_string, *) "No Gamma file found"
      LOG_ERROR(log_string)
    END IF
  END IF

END SUBROUTINE GET_SAFE_GAMMA_SC

SUBROUTINE GET_GAMMA_SC(Gamma_SC, path)
  !! This subroutine reads the Gamma_SC file from self-consistent simulation at given path.
  !! User should check if the file exists before calling this subroutine
  CHARACTER(LEN=*), INTENT(IN) :: path !! path from which to read Gamma file
  COMPLEX*16, INTENT(OUT) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS) !! Gamma to be filled
  INTEGER*4 :: n, lat, orb, spin, band
  INTEGER*4 :: n_read, lat_read, orb_read, spin_read, band_read
  REAL*8 :: Gamma_re, Gamma_im
  CHARACTER(LEN=20) :: output_format

#ifndef READ_OLD
  output_format = '(5I5, 2E15.5)'
#else
  output_format = '(4I5, 2E15.5)'
#endif

  OPEN (unit=9, FILE=path, FORM="FORMATTED", ACTION="READ", STATUS="OLD")
  READ (9, *)
  DO band = 1, SUBBANDS
    DO spin = 1, 2
      DO n = 1, N_NEIGHBOURS
        DO lat = 1, LAYER_COUPLINGS
          DO orb = 1, ORBITALS
#ifndef READ_OLD
            READ (9, output_format) band_read, spin_read, n_read, lat_read, orb_read, Gamma_re, Gamma_im
            Gamma_SC(orb_read, n_read, spin_read, lat_read, band_read) = DCMPLX(Gamma_re, Gamma_im) * meV2au
#else
            READ (9, output_format) spin_read, n_read, lat_read, orb_read, Gamma_re, Gamma_im
            Gamma_SC(orb_read, n_read, spin_read, lat_read, band) = DCMPLX(Gamma_re, Gamma_im) * meV2au
#endif
          END DO
        END DO
        READ (9, *)
        READ (9, *)
      END DO
      DO n = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
        DO lat = 1, SUBLATTICES
          DO orb = 1, ORBITALS
#ifndef READ_OLD
            READ (9, output_format) band_read, spin_read, n_read, lat_read, orb_read, Gamma_re, Gamma_im
            Gamma_SC(orb_read, n_read, spin_read, lat_read, band_read) = DCMPLX(Gamma_re, Gamma_im) * meV2au
#else
            READ (9, output_format) spin_read, n_read, lat_read, orb_read, Gamma_re, Gamma_im
            Gamma_SC(orb_read, n_read, spin_read, lat_read, band) = DCMPLX(Gamma_re, Gamma_im) * meV2au
#endif
          END DO
        END DO
        READ (9, *)
        READ (9, *)
      END DO
    END DO
  END DO
  CLOSE (9)

END SUBROUTINE GET_GAMMA_SC

SUBROUTINE GET_SAFE_CHARGE_DENS(Charge_dens, input_path)
  !! This subroutine tries to read the Chargen_dens file from self-consistent simulation, checking for final result first,
  !! then checking the latest iteration file, and logging an error if neither is found.
  REAL*8, INTENT(OUT) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS) !! Charge density to be filled
  CHARACTER(LEN=*), INTENT(IN) :: input_path !! Input path in which to look for Charge file
  LOGICAL :: file_exists !! Whether the file exists or not
  !Read charges
  INQUIRE (FILE=TRIM(input_path)//"OutputData/Charge_dens_final.dat", EXIST=file_exists)
  IF (file_exists) THEN
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(input_path)//"OutputData/Charge_dens_final.dat")
  ELSE
    INQUIRE (FILE=TRIM(input_path)//"OutputData/Charge_dens_iter.dat", EXIST=file_exists)
    IF (file_exists) THEN
      CALL GET_CHARGE_DENS(Charge_dens, TRIM(input_path)//"OutputData/Charge_dens_iter.dat")
    ELSE
      WRITE (log_string, *) "No charge density file found"
      LOG_ERROR(log_string)
    END IF
  END IF
END SUBROUTINE GET_SAFE_CHARGE_DENS

SUBROUTINE GET_CHARGE_DENS(Charge_dens, path)
  !! Read a file containing the charge density from given path. Whether file exists has to be checked before.
  CHARACTER(LEN=*), INTENT(IN) :: path !! Path to file
  REAL*8, INTENT(OUT) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS) !! Charge density to be filled
  INTEGER*4 :: spin, lat, orb, n, band, band_read
  CHARACTER(LEN=20) :: output_format

#ifndef READ_OLD
  output_format = '(4I5, 1E15.5)'
#else
  output_format = '(3I5, 1E15.5)'
#endif

  OPEN (unit=9, FILE=path, FORM="FORMATTED", ACTION="READ", STATUS="OLD")
  READ (9, *)
  DO band = 1, SUBBANDS
    DO n = 1, DIM_POSITIVE_K
#ifndef READ_OLD
      READ (9, output_format) band_read, spin, lat, orb, Charge_dens(n, band)
#else
      READ (9, output_format) spin, lat, orb, Charge_dens(n, band)
#endif
    END DO
  END DO

  CLOSE (9)

END SUBROUTINE GET_CHARGE_DENS

SUBROUTINE SET_HAMILTONIAN_PARAMS()
    !! This subroutine sets global variables that define dimension of hamiltonian in calculation
    !! SUBLATTICES should have been set beffore at reading input.nml
  ORBITALS = 3
  TBA_DIM = ORBITALS * SUBLATTICES
  DIM_POSITIVE_K = TBA_DIM * 2  !Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
  DIM = DIM_POSITIVE_K * 2    !*2 to transform to Nambu Space.
  LAYER_COUPLINGS = 2 * (SUBLATTICES - 1)
END SUBROUTINE

END MODULE mod_reader
