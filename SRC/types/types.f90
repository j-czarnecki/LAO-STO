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

TYPE crs_matrix_t
  !! Compressed Row Storage matrix type.
  INTEGER(INT32) :: n_nonzero !! Number of nonzero elements in the matrix.
  REAL(REAL64), ALLOCATABLE :: Values(:) !! Nonzero values in the matrix.
  INTEGER(INT32), ALLOCATABLE :: Column_indices(:) !! Column indeces of consecutive nonzero values.
  INTEGER(INT32), ALLOCATABLE :: Row_indices(:) !! Indeces in Values/Column_indeces arrays that start a new row.
END TYPE crs_matrix_t

TYPE derived_dimensions_t
  !! Parameters derived from degrees of freedom and discretization specification.
  REAL(REAL64) :: dr_k = 0.0d0 !! Radial step in k-space for integration.
  REAL(REAL64) :: dphi_k = 0.0d0 !! Angular step in k-space for integration.
  INTEGER(INT32) :: TBA_DIM = 0 !! Tight-binding Hamiltonian dimension (without spin and Nambu space).
  INTEGER(INT32) :: DIM_POSITIVE_K = 0  !! Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
  INTEGER(INT32) :: DIM = 0    !! *2 to transform to Nambu Space.
  INTEGER(INT32) :: LAYER_COUPLINGS = 0 !! This determines how many layer-related superconducting parameters have to be calculated.
END TYPE derived_dimensions_t

TYPE discretization_t
  !! Discretization parameters for k-space integration and degrees of freedom.
  !! Later on used to set derived_dimensions_t parameters.
  INTEGER(INT32) :: k1_steps = 0
  INTEGER(INT32) :: k2_steps = 0
  INTEGER(INT32) :: SUBLATTICES = 2
  INTEGER(INT32) :: SUBBANDS = 1
  INTEGER(INT32) :: ORBITALS = 3
  !Derived
  TYPE(derived_dimensions_t) :: derived
END TYPE discretization_t

TYPE subband_params_t
  !! Physical parameters for a single subband.
  REAL(REAL64) :: t_D = 0.0d0 !! Direct hopping (sigma overlap) energy betweend t_{2g} d orbitals [meV].
  REAL(REAL64) :: t_I = 0.0d0 !! Indirect hopping (pi overlap) energy betweend t_{2g} d orbitals [meV].
  REAL(REAL64) :: t_Rashba = 0.0d0 !! Rashba spin-orbit coupling strength [meV].
  REAL(REAL64) :: lambda_SOC = 0.0d0 !! Atomic spin-orbit coupling strength [meV].
  REAL(REAL64) :: delta_trigonal = 0.0d0 !! Trigonal crystal field splitting [meV].
  REAL(REAL64) :: zeta_tetragonal = 0.0d0 !! Tetragonal crystal field splitting [meV].
  INTEGER(INT32) :: orb_affected_tetragonal = 1 !! Orbital affected by tetragonal distortion: 1 - d_yz, 2 - d_zx, 3 - d_xy.
  REAL(REAL64) :: v = 0.0d0 !! Electric potential gradient at the interface [meV per sublattice distance].
  REAL(REAL64) :: V_pdp = 0.0d0 !! Slater-Koster integral for p-d orbitals hybrdization via pi overlap [meV].
  REAL(REAL64) :: V_pds = 0.0d0 !! Slater-Koster integral for p-d orbitals hybrdization via sigma overlap [meV].
  TYPE(crs_matrix_t) :: J_tensor !! Electron-hole coupling tensor in CRS format.
  REAL(REAL64) :: U_HUB = 0.0d0 !! Hubbard on-site intraorbital repulsion energy [meV].
  REAL(REAL64) :: V_HUB = 0.0d0 !! Hubbard on-site interorbital repulsion energy [meV].
  REAL(REAL64) :: E_Fermi = 0.0d0 !! Fermi energy [meV].
  REAL(REAL64), ALLOCATABLE :: V_layer(:) !! Layer-dependent potential energy [meV].
  REAL(REAL64), ALLOCATABLE :: Subband_energies(:) !! Subband bottom energies [meV].
  REAL(REAL64) :: g_factor = 0.0d0 !! Effective g-factor for Zeeman coupling.
  !Derived
  REAL(REAL64) :: eta_p = 0.0d0 !! Parameter defining strength of p-d hybridization.
END TYPE subband_params_t

TYPE external_params_t
  !! External physical parameters.
  REAL(REAL64) :: T = 0.0d0 !! Temparature
  REAL(REAL64) :: B_field(3) !! External magnetic field vector [T].
END TYPE external_params_t

TYPE physical_params_t
  TYPE(subband_params_t) :: subband_params
  TYPE(external_params_t) :: external
END TYPE physical_params_t

TYPE self_consistency_t
  !! Parameters for self-consistent calculation of superconducting order parameter and charge.
  LOGICAL :: read_gamma_from_file = .FALSE. !! Whether to read initial gamma from file.
  CHARACTER(1000) :: path_to_gamma_start !! Path to file with initial gamma values. Has to end with '/'.
  LOGICAL :: read_charge_from_file = .FALSE. !! Whether to read initial charge from file.
  CHARACTER(1000) :: path_to_charge_start !! Path to file with initial charge values. Has to end with '/'.
  REAL(REAL64) :: gamma_start = 0.0d0 !! Initial value of gamma for all layer couplings if not read from file.
  REAL(REAL64) :: gamma_nnn_start = 0.0d0 !! Initial value of gamma for next-nearest-neighbor layer couplings if not read from file.
  REAL(REAL64) :: charge_start = 0.0d0 !! Initial value of charge for all layers if not read from file.
  INTEGER(INT32) :: max_sc_iter = 0 !! Maximum number of self-consistency iterations.
  REAL(REAL64) :: sc_alpha = 0.0d0 !! Mixing parameter for self-consistency iterations.
  REAL(REAL64) :: sc_alpha_adapt = 0.0d0 !! Damping factor for mixing parameter in self-consistency iterations.
                                         !! Used when subsequent iterations are divergent.
  REAL(REAL64) :: gamma_eps_convergence = 0.0d0 !! Absolute error for convergence of superconducting order parameter.
  REAL(REAL64) :: charge_eps_convergence = 0.0d0 !! Absolute error for convergence of charge.
END TYPE self_consistency_t

TYPE romberg_integration_t
  !! Parameters for Romberg integration over the Brillouin zone.
  REAL(REAL64) :: romb_eps_x = 0.0d0 !! Desired relative accuracy for Romberg integration in k1 direction.
  INTEGER(INT32) :: interpolation_deg_x = 0 !! Degree of polynomial interpolation in k1 direction.
  INTEGER(INT32) :: max_grid_refinements_x = 0 !! Maximum number of grid refinements in k1 direction.
                                               !! Each refinement doubles the number of grid points.
  REAL(REAL64) :: romb_eps_y = 0.0d0 !! Desired relative accuracy for Romberg integration in k2 direction.
  INTEGER(INT32) :: interpolation_deg_y = 0 !! Degree of polynomial interpolation in k2 direction.
  INTEGER(INT32) :: max_grid_refinements_y = 0 !! Maximum number of grid refinements in k2 direction.
                                               !! Each refinement doubles the number of grid points.
END TYPE romberg_integration_t

TYPE sc_input_params_t
  !! Input parameters for self-consistent superconductivity calculations.
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
  LOGICAL :: enable = .FALSE. !! Whether to perform SC gap calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  REAL(REAL64) :: dE = 0.0d0 !! Energy window around Fermi level for SC gap calculation [meV].
  INTEGER(INT32) :: Nk_points = 0 !! Number of k-points along each direction in the Brillouin zone for initial
                                  !! search for SC gap states.
  INTEGER(INT32) :: Nk_points_refined = 0 !! Number of k-points along each direction in the Brillouin zone
                                          !! for refined SC gap calculation.
                                          !! Whenever initial search finds states within dE of E_Fermi,
                                          !! refined grid is used to calculate SC gap more accurately
                                          !! around this point.
END TYPE post_sc_gap_t

TYPE post_chern_number_t
  !! Chern number calculation
  LOGICAL :: enable = .FALSE. !! Whether to perform Chern number calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  INTEGER(INT32) :: Nk_points = 0 !! Number of k-points along each direction in the Brillouin
                                  !! zone for Chern number calculation.
END TYPE post_chern_number_t

TYPE post_dispersion_relation_t
  !! Dispersion relation calculation
  LOGICAL :: enable = .FALSE. !! Whether to perform dispersion relation calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  LOGICAL :: include_sc = .FALSE. !! Whether to include hole bands in dispersion relation.
  INTEGER(INT32) :: Nr_points = 0 !! Number of radial k-points along each direction in the Brillouin zone.
  INTEGER(INT32) :: Nphi_points = 0 !! Number of angular k-points in the Brillouin zone.
END TYPE post_dispersion_relation_t

TYPE post_dos_t
  !! DOS calculation
  LOGICAL :: enable = .FALSE. !! Whether to perform DOS calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  REAL(REAL64) :: E_min = 0 !! Minimum energy for DOS calculation [meV].
  REAL(REAL64) :: E_max = 0 !! Maximum energy for DOS calculation [meV].
  REAL(REAL64) :: dE0 = 0 !! Initial energy resolution for DOS calculation [meV].
  REAL(REAL64) :: zeta_DOS = 0 !! Broadening parameter for Dirac's delta function for DOS calculation [meV].
  LOGICAL :: include_sc = .FALSE. !! Whether to include hole bands in DOS calculation.
  INTEGER(INT32) :: Nk_points = 0 !! Number of k-points along each direction in the Brillouin zone for initial
                                  !! DOS calculation.
  INTEGER(INT32) :: Nk_points_refined = 0 !! Number of k-points along each direction in the Brillouin zone
                                          !! for refined DOS calculation.
                                          !! Whenever initial DOS calculation finds states within zeta_DOS
                                          !! of E_Fermi, refined grid is used to calculate DOS more accurately
                                          !! around this point.
END TYPE post_dos_t

TYPE post_gamma_k_t
  !! Calculations o momentum dependent order parameter gamma(k).
  LOGICAL :: enable = .FALSE. !! Whether to perform gamma(k) calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  INTEGER(INT32) :: Nk_points = 0 !! Number of k-points along each direction in the Brillouin zone for gamma(k) calculation.
END TYPE post_gamma_k_t

TYPE post_projections_t
  !! Projections calculation using integration over Brillouin zone.
  LOGICAL :: enable = .FALSE. !! Whether to perform projections calculation.
  CHARACTER(1000) :: path = "" !! Path from which to read Gamma and Charge.
  INTEGER(INT32) :: Nr_points = 0 !! Number of radial k-points along each direction in the Brillouin zone.
  INTEGER(INT32) :: Nphi_points = 0 !! Number of angular k-points in the Brillouin zone.
END TYPE post_projections_t

TYPE post_input_params_t
  !! Input parameters for postprocessing calculations.
  TYPE(post_sc_gap_t) :: sc_gap
  TYPE(post_chern_number_t) :: chern
  TYPE(post_dispersion_relation_t) :: dispersion
  TYPE(post_dos_t) :: dos
  TYPE(post_gamma_k_t) :: gamma_k
  TYPE(post_projections_t) :: projections
END TYPE post_input_params_t

END MODULE types
