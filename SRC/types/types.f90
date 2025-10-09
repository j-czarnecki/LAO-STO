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

MODULE types
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
IMPLICIT NONE

TYPE derived_dimensions_t
  !Those parameters are in fact derived and recalculated if needed at GET_INPUT
  !(see SET_HAMILTONIAN_PARAMS)
  REAL(REAL64) :: dr_k = 0.0d0
  REAL(REAL64) :: dphi_k = 0.0d0
  INTEGER(INT32) :: TBA_DIM = 0
  INTEGER(INT32) :: DIM_POSITIVE_K = 0  !Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
  INTEGER(INT32) :: DIM = 0    !*2 to transform to Nambu Space.
  INTEGER(INT32) :: LAYER_COUPLINGS = 0 !This determines how many layer-related superconducting parameters have to be calculated.

END TYPE derived_dimensions_t

TYPE discretization_t
  !Discretization
  INTEGER(INT32) :: k1_steps = 0
  INTEGER(INT32) :: k2_steps = 0
  INTEGER(INT32) :: SUBLATTICES = 2
  INTEGER(INT32) :: SUBBANDS = 1
  INTEGER(INT32) :: ORBITALS = 3
  !Derived
  TYPE(derived_dimensions_t) :: derived
END TYPE discretization_t

TYPE subband_params_t
  REAL(REAL64) :: t_D = 0.0d0
  REAL(REAL64) :: t_I = 0.0d0
  REAL(REAL64) :: t_Rashba = 0.0d0
  REAL(REAL64) :: lambda_SOC = 0.0d0
  REAL(REAL64) :: delta_trigonal = 0.0d0
  REAL(REAL64) :: zeta_tetragonal = 0.0d0
  INTEGER(INT32) :: orb_affected_tetragonal = 1
  REAL(REAL64) :: v = 0.0d0
  REAL(REAL64) :: V_pdp = 0.0d0
  REAL(REAL64) :: V_pds = 0.0d0
  REAL(REAL64) :: J_SC_tensor(SPINS, SPINS, SPINS, SPINS) = 0.0d0
  REAL(REAL64) :: nearest_interorb_multiplier = 0.0d0
  REAL(REAL64) :: J_SC_NNN_tensor(SPINS, SPINS, SPINS, SPINS) = 0.0d0
  REAL(REAL64) :: next_interorb_multiplier = 0.0d0
  REAL(REAL64) :: U_HUB = 0.0d0
  REAL(REAL64) :: V_HUB = 0.0d0
  REAL(REAL64) :: E_Fermi = 0.0d0
  REAL(REAL64), ALLOCATABLE :: V_layer(:)
  REAL(REAL64), ALLOCATABLE :: Subband_energies(:)
  REAL(REAL64) :: g_factor = 0.0d0
  !Derived
  REAL(REAL64) :: eta_p = 0.0d0
END TYPE subband_params_t

TYPE external_params_t
  REAL(REAL64) :: T = 0.0d0 !! Temparature
  REAL(REAL64) :: B_field(3)
END TYPE external_params_t

TYPE physical_params_t
  TYPE(subband_params_t) :: subband_params
  TYPE(external_params_t) :: external
END TYPE physical_params_t

TYPE self_consistency_t
  LOGICAL :: read_gamma_from_file = .FALSE.
  CHARACTER(1000) :: path_to_gamma_start
  LOGICAL :: read_charge_from_file = .FALSE.
  CHARACTER(1000) :: path_to_charge_start
  REAL(REAL64) :: gamma_start = 0.0d0
  REAL(REAL64) :: gamma_nnn_start = 0.0d0
  REAL(REAL64) :: charge_start = 0.0d0
  INTEGER(INT32) :: max_sc_iter = 0
  REAL(REAL64) :: sc_alpha = 0.0d0
  REAL(REAL64) :: sc_alpha_adapt = 0.0d0
  REAL(REAL64) :: gamma_eps_convergence = 0.0d0
  REAL(REAL64) :: charge_eps_convergence = 0.0d0
END TYPE self_consistency_t

TYPE romberg_integration_t
  REAL(REAL64) :: romb_eps_x = 0.0d0
  INTEGER(INT32) :: interpolation_deg_x = 0
  INTEGER(INT32) :: max_grid_refinements_x = 0
  REAL(REAL64) :: romb_eps_y = 0.0d0
  INTEGER(INT32) :: interpolation_deg_y = 0
  INTEGER(INT32) :: max_grid_refinements_y = 0
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
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  REAL(REAL64) :: dE = 0.0d0
  INTEGER(INT32) :: Nk_points = 0
  INTEGER(INT32) :: Nk_points_refined = 0
END TYPE post_sc_gap_t

TYPE post_chern_number_t
    !! Chern number calculation
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  INTEGER(INT32) :: Nk_points = 0
END TYPE post_chern_number_t

TYPE post_dispersion_relation_t
    !! Dispersion relation calculation
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  LOGICAL :: include_sc = .FALSE.
  INTEGER(INT32) :: Nr_points = 0
  INTEGER(INT32) :: Nphi_points = 0
END TYPE post_dispersion_relation_t

TYPE post_dos_t
  !! DOS calculation
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  REAL(REAL64) :: E_min = 0
  REAL(REAL64) :: E_max = 0
  REAL(REAL64) :: dE0 = 0
  REAL(REAL64) :: zeta_DOS = 0
  LOGICAL :: include_sc = .FALSE.
  INTEGER(INT32) :: Nk_points = 0
  INTEGER(INT32) :: Nk_points_refined = 0
END TYPE post_dos_t

TYPE post_gamma_k_t
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  INTEGER(INT32) :: Nk_points = 0
END TYPE post_gamma_k_t

TYPE post_projections_t
  LOGICAL :: enable = .FALSE.
  CHARACTER(1000) :: path = ""
  INTEGER(INT32) :: Nr_points = 0
  INTEGER(INT32) :: Nphi_points = 0
END TYPE post_projections_t

TYPE post_input_params_t
  TYPE(post_sc_gap_t) :: sc_gap
  TYPE(post_chern_number_t) :: chern
  TYPE(post_dispersion_relation_t) :: dispersion
  TYPE(post_dos_t) :: dos
  TYPE(post_gamma_k_t) :: gamma_k
  TYPE(post_projections_t) :: projections
END TYPE post_input_params_t

END MODULE types
