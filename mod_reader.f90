MODULE mod_reader
USE mod_parameters
IMPLICIT NONE
SAVE

!Default values, overwritten in get_input

!Physical parameters
REAL*8 :: T = 0.
REAL*8 :: t_D = 0.
REAL*8 :: t_I = 0.
REAL*8 :: t_Rashba = 0.
REAL*8 :: lambda_SOC = 0.
REAL*8 :: DELTA_TRI = 0.
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

!Discretization
INTEGER*4 :: k1_steps = 0
INTEGER*4 :: k2_steps = 0



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
REAL*8 eta_p
REAL*8 :: dk1, dk2, domega


!Used for postprocessing
!Superconducting gap calculation
LOGICAL :: enable_sc_gap_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_sc_gap = ""
REAL*8 :: dE_sc_gap = 0.
INTEGER*4 :: Nk_points_sc_gap = 0

!For chern number calculation
LOGICAL :: enable_chern_number_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_chern_number = ""
INTEGER*4 :: Nk_points_chern_number = 0

!Dispersion relation calculation
LOGICAL :: enable_dispersion_relation_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_dispersion_relation = ""
INTEGER*4 :: Nk_points_dispersion_relation = 0


!DOS calculation
LOGICAL :: enable_dos_calc = .FALSE.
CHARACTER(1000) :: path_to_run_dir_dos = ""
REAL*8 :: E_DOS_min = 0
REAL*8 :: E_DOS_max = 0
REAL*8 :: dE0 = 0
REAL*8 :: zeta_DOS = 0
INTEGER*4 :: Nk_points_dos = 0


NAMELIST /physical_params/  &
& T,                        &
& t_D,                      &
& t_I,                      &
& t_Rashba,                 &
& lambda_SOC,               &
& DELTA_TRI,                &
& v,                        &
& V_pdp,                    &
& V_pds,                    &
& J_SC,                     &
& J_SC_PRIME,               &
& J_SC_NNN,                 &
& J_SC_PRIME_NNN,           &
& U_HUB,                    &
& V_HUB,                    &
& E_Fermi

NAMELIST /discretization/ &
& k1_steps,               &
& k2_steps

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
& Nk_points_sc_gap

NAMELIST /chern_number_calculation/ &
& enable_chern_number_calc, &
& path_to_run_dir_chern_number, &
& Nk_points_chern_number

NAMELIST /dispersion_relation_calculation/ &
& enable_dispersion_relation_calc, &
& path_to_run_dir_dispersion_relation, &
& Nk_points_dispersion_relation

NAMELIST /dos_calculation/ &
& enable_dos_calc, &
& path_to_run_dir_dos, &
& E_DOS_min, &
& E_DOS_max, &
& dE0, &
& zeta_DOS, &
& Nk_points_dos

CONTAINS
SUBROUTINE GET_INPUT(nmlfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: nmlfile

    OPEN(unit = 9, FILE=nmlfile, FORM = "FORMATTED", ACTION = "READ", STATUS="OLD")

    !TODO: WRITE BETTER CHECKS!!!!!!!!!!!!!!!!!!
    READ(9,NML=physical_params)
    !Check input data
    IF (T < 0) STOP "Temperature in kelvins must be >= 0!"

    !Change to atomic units
    t_D = t_D * meV2au
    t_I = t_I * meV2au
    t_Rashba = t_Rashba * meV2au
    lambda_SOC = lambda_SOC * meV2au
    DELTA_TRI = DELTA_TRI * meV2au
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

    READ(9,NML=discretization)
    IF ((k1_steps .LE. 0) .OR. (k2_steps .LE. 0)) STOP "k_steps must be > 0"

    READ(9,NML=self_consistency)
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
    dk1 = K1_MAX / k1_steps
    dk2 = K2_MAX / k2_steps
    domega = ABS(dk1*dk2*SIN(2*PI/3.))/(SIN(2*PI/3.)*K1_MAX*K2_MAX)
    eta_p = v * SQRT(3.) / 3.905 * nm2au

    READ(9, NML=romberg_integration)
    IF (romb_eps_x .LE. 0) STOP "romb_eps_x must be > 0"
    IF (interpolation_deg_x .LE. 0) STOP "interpolation_deg_x must be > 0"
    IF (max_grid_refinements_x .LE. 0) STOP "max_grid_refinements_x must be > 0"

    IF (romb_eps_y .LE. 0) STOP "romb_eps_y must be > 0"
    IF (interpolation_deg_y .LE. 0) STOP "interpolation_deg_y must be > 0"
    IF (max_grid_refinements_y .LE. 0) STOP "max_grid_refinements_y must be > 0"



    CLOSE(9)

END SUBROUTINE GET_INPUT

SUBROUTINE GET_POSTPROCESSING_INPUT(nmlfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: nmlfile

    OPEN(unit = 9, FILE=nmlfile, FORM = "FORMATTED", ACTION = "READ", STATUS="OLD")
    READ(9, NML = sc_gap_calculation)
    IF (enable_sc_gap_calc) THEN
        IF (Nk_points_sc_gap .LE. 0) STOP "Nk_points_sc_gap must be > 0"
        IF (path_to_run_dir_sc_gap == "") STOP "path_to_run_dir_sc_gap must not be empty"
        IF (dE_sc_gap .LE. 0) STOP "dE_sc_gap must be > 0"

        dE_sc_gap = dE_sc_gap * meV2au
    END IF

    READ(9, NML = chern_number_calculation)
    IF (enable_chern_number_calc) THEN
        IF (Nk_points_chern_number .LE. 0) STOP "Nk_points_chern_number must be > 0"
        IF (path_to_run_dir_chern_number == "") STOP "path_to_run_dir_chern_number must not be empty"
    END IF

    READ(9, NML = dispersion_relation_calculation)
    IF (enable_dispersion_relation_calc) THEN
        IF (Nk_points_dispersion_relation .LE. 0) STOP "Nk_points_dispersion_relation must be > 0"
        IF (path_to_run_dir_dispersion_relation == "") STOP "path_to_run_dir_dispersion_relation must not be empty"
    END IF

    READ(9, NML = dos_calculation)
    IF (enable_dos_calc) THEN
        IF (path_to_run_dir_dos == "") STOP "path_to_run_dir_dos must not be empty"
        IF (E_DOS_min .GT. E_DOS_max) STOP "E_DOS_min must be smaller than E_DOS_max"
        IF (dE0 .LE. 0) STOP "dE0 must be > 0"
        IF (zeta_DOS .LE. 0) STOP "zeta_DOS must be > 0"
        IF (Nk_points_dos .LE. 0) STOP "Nk_points_dos must be > 0"
        E_DOS_min = E_DOS_min * meV2au
        E_DOS_max = E_DOS_max * meV2au
        dE0 = dE0 * meV2au
    END IF

    CLOSE(9)

END SUBROUTINE GET_POSTPROCESSING_INPUT


SUBROUTINE GET_GAMMA_SC(Gamma_SC, path)
    CHARACTER(LEN=*), INTENT(IN) :: path
    COMPLEX*16, INTENT(OUT) :: Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES)
    INTEGER*4 :: n, lat, orb,spin
    INTEGER*4 :: n_read, lat_read, orb_read, spin_read
    REAL*8 :: Gamma_re, Gamma_im
    CHARACTER(LEN=20) :: output_format

    output_format = '(4I5, 2E15.5)'

    OPEN(unit = 9, FILE=path, FORM = "FORMATTED", ACTION = "READ", STATUS="OLD")
    READ(9,*)

    DO spin =1, 2
        DO n = 1, N_ALL_NEIGHBOURS
            DO lat = 1, SUBLATTICES
                DO orb = 1, ORBITALS
                    READ(9, output_format) spin_read, n_read, lat_read, orb_read, Gamma_re, Gamma_im
                    Gamma_SC(orb_read, n_read, spin_read, lat_read) = DCMPLX(Gamma_re, Gamma_im)*meV2au
                END DO
            END DO
            READ(9,*)
            READ(9,*)
        END DO
    END DO
    CLOSE(9)

END SUBROUTINE GET_GAMMA_SC

SUBROUTINE GET_CHARGE_DENS(Charge_dens, path)
    CHARACTER(LEN=*), INTENT(IN) :: path
    REAL*8, INTENT(OUT) :: Charge_dens(DIM_POSITIVE_K)
    INTEGER*4 :: spin, lat, orb, n
    CHARACTER(LEN=20) :: output_format

    output_format = '(3I5, 1E15.5)'

    OPEN(unit = 9, FILE=path, FORM = "FORMATTED", ACTION = "READ", STATUS="OLD")
    READ(9,*)

    DO n = 1, DIM_POSITIVE_K
        READ(9, output_format) spin, lat, orb, Charge_dens(n)
    END DO

    CLOSE(9)

END SUBROUTINE GET_CHARGE_DENS

END MODULE mod_reader
