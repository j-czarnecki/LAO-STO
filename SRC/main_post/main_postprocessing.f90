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
PROGRAM main_postprocessing
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE hamiltonians
USE parameters
USE utilities
USE writers
USE reader
USE local_integrand
USE postprocessing
USE logger
USE types
IMPLICIT NONE

TYPE(post_input_params_t) :: post_input_params

CALL INIT_LOGGER("postprocessing")

WRITE (log_string, *) "Reading postprocessing input"
LOG_INFO(log_string)

CALL GET_POSTPROCESSING_INPUT('./postprocessing_input.nml', post_input_params)

IF (post_input_params % dispersion % enable) THEN
  WRITE (log_string, *) "Calculating dispersion relation"
  LOG_INFO(log_string)
  CALL CALCULATE_DISPERSION(post_input_params % dispersion)
END IF

IF (post_input_params % dos % enable) THEN
  WRITE (log_string, *) "Calculating DOS"
  LOG_INFO(log_string)
  CALL CALCULATE_DOS(post_input_params % dos)
END IF

IF (post_input_params % chern % enable) THEN
  WRITE (log_string, *) "Calculating Chern number"
  LOG_INFO(log_string)
  CALL CALCULATE_CHERN_PARAMS(post_input_params % chern)
END IF

IF (post_input_params % sc_gap % enable) THEN
  WRITE (log_string, *) "Calculating superconducting gap"
  LOG_INFO(log_string)
  CALL CALCULATE_SUPERCONDUCTING_GAP(post_input_params % sc_gap)
END IF

IF (post_input_params % gamma_k % enable) THEN
  WRITE (log_string, *) "Calculating Gamma_K"
  LOG_INFO(log_string)
  CALL CALCULATE_GAMMA_K(post_input_params % gamma_k)
END IF

IF (post_input_params % projections % enable) THEN
  WRITE (log_string, *) "Calculating projections"
  LOG_INFO(log_string)
  CALL CALCULATE_PROJECTIONS(post_input_params % projections)
END IF

CALL CLOSE_LOGGER()

END PROGRAM main_postprocessing
