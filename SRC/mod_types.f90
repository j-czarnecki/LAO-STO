MODULE mod_types

TYPE derived_dimensions_t
  !Those parameters are in fact derived and recalculated if needed at GET_INPUT
  !(see SET_HAMILTONIAN_PARAMS)
  INTEGER*4 :: ORBITALS = 0
  INTEGER*4 :: TBA_DIM = 0
  INTEGER*4 :: DIM_POSITIVE_K = 0  !Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
  INTEGER*4 :: DIM = 0    !*2 to transform to Nambu Space.
  INTEGER*4 :: LAYER_COUPLINGS = 0 !This determines how many layer-related superconducting parameters have to be calculated.

END TYPE derived_dimensions_t

TYPE discretization_t
  !Discretization
  INTEGER*4 :: k1_steps = 0
  INTEGER*4 :: k2_steps = 0
  INTEGER*4 :: SUBLATTICES = 2
  INTEGER*4 :: SUBBANDS = 1
  !Derived
  REAL*8 :: dk1 = 0.0d0
  REAL*8 :: dk2 = 0.0d0
  REAL*8 :: domega = 0.0d0
  TYPE(derived_dimensions_t) :: derived
END TYPE discretization_t

TYPE subband_params_t
  REAL*8 :: t_D = 0.0d0
  REAL*8 :: t_I = 0.0d0
  REAL*8 :: t_Rashba = 0.0d0
  REAL*8 :: lambda_SOC = 0.0d0
  REAL*8 :: DELTA_TRI = 0.0d0
  REAL*8 :: v = 0.0d0
  REAL*8 :: V_pdp = 0.0d0
  REAL*8 :: V_pds = 0.0d0
  REAL*8 :: J_SC = 0.0d0
  REAL*8 :: J_SC_PRIME = 0.0d0
  REAL*8 :: J_SC_NNN = 0.0d0
  REAL*8 :: J_SC_PRIME_NNN = 0.0d0
  REAL*8 :: U_HUB = 0.0d0
  REAL*8 :: V_HUB = 0.0d0
  REAL*8 :: E_Fermi = 0.0d0
  REAL*8 :: V_layer = 0.0d0
  REAL*8 :: Subband_energies = 0.0d0
  !Derived
  REAL*8 :: eta_p = 0.0d0
END TYPE subband_params_t

TYPE external_params_t
  REAL*8 :: B_field(3) = (/0.0d0, 0.0d0, 0.0d0/)
END TYPE external_params_t

TYPE physical_params_t
  REAL*8 :: T = 0.0d0                         !! Temparature
  TYPE(subband_params_t) :: subband_params
  TYPE(external_params_t) :: external
END TYPE physical_params_t

TYPE self_consistency_t
  LOGICAL :: read_gamma_from_file = .FALSE.
  CHARACTER(1000) :: path_to_gamma_start
  LOGICAL :: read_charge_from_file = .FALSE.
  CHARACTER(1000) :: path_to_charge_start
  REAL*8 :: gamma_start = 0.0d0
  REAL*8 :: gamma_nnn_start = 0.0d0
  REAL*8 :: charge_start = 0.0d0
  INTEGER*4 :: max_sc_iter = 0
  REAL*8 :: sc_alpha = 0.0d0
  REAL*8 :: sc_alpha_adapt = 0.0d0
  REAL*8 :: gamma_eps_convergence = 0.0d0
  REAL*8 :: charge_eps_convergence = 0.0d0
END TYPE self_consistency_t

TYPE romberg_integration_t
  REAL*8 :: romb_eps_x = 0.0d0
  INTEGER*4 :: interpolation_deg_x = 0
  INTEGER*4 :: max_grid_refinements_x = 0
  REAL*8 :: romb_eps_y = 0.0d0
  INTEGER*4 :: interpolation_deg_y = 0
  INTEGER*4 :: max_grid_refinements_y = 0
END TYPE romberg_integration_t

TYPE sc_input_params_t
  TYPE(discretization_t) :: discretization
  TYPE(physical_params_t) :: physical
  TYPE(self_consistency_t) :: self_consistency
  TYPE(romberg_integration_t) :: romberg
END TYPE sc_input_params_t

! ----------------------------------------------------------------------
! ---------------------- Used for postprocessing -----------------------
! ----------------------------------------------------------------------
TYPE post_sc_gap_t
    !! Superconducting gap calculation
  LOGICAL :: enable_sc_gap_calc = .FALSE.
  CHARACTER(1000) :: path_to_run_dir_sc_gap = ""
  REAL*8 :: dE_sc_gap = 0.0d0
  INTEGER*4 :: Nk_points_sc_gap = 0
  INTEGER*4 :: Nk_points_sc_gap_refined = 0
END TYPE post_sc_gap_t

TYPE post_chern_number_t
    !! Chern number calculation
  LOGICAL :: enable_chern_number_calc = .FALSE.
  CHARACTER(1000) :: path_to_run_dir_chern_number = ""
  INTEGER*4 :: Nk_points_chern_number = 0
END TYPE post_chern_number_t

TYPE post_dispersion_relation_t
    !! Dispersion relation calculation
  LOGICAL :: enable_dispersion_relation_calc = .FALSE.
  CHARACTER(1000) :: path_to_run_dir_dispersion_relation = ""
  LOGICAL :: include_sc_in_dispersion = .FALSE.
  INTEGER*4 :: Nk_points_dispersion_relation = 0
END TYPE post_dispersion_relation_t

TYPE post_dos_t
    !! DOS calculation
  LOGICAL :: enable_dos_calc = .FALSE.
  CHARACTER(1000) :: path_to_run_dir_dos = ""
  REAL*8 :: E_DOS_min = 0
  REAL*8 :: E_DOS_max = 0
  REAL*8 :: dE0 = 0
  REAL*8 :: zeta_DOS = 0
  LOGICAL :: include_sc_in_dos = .FALSE.
  INTEGER*4 :: Nk_points_dos = 0
  INTEGER*4 :: Nk_points_dos_refined = 0
END TYPE post_dos_t

TYPE post_input_params_t
  TYPE(post_sc_gap_t) :: post_sc
  TYPE(post_chern_number_t) :: post_chern
  TYPE(post_dispersion_relation_t) :: post_dispersion
  TYPE(post_dos_t) :: post_dos
END TYPE post_input_params_t

END MODULE mod_types
