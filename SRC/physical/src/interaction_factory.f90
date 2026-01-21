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

MODULE interaction_factory
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE utilities
USE parameters
USE types
USE logger
IMPLICIT NONE
CONTAINS

SUBROUTINE ASSIGN_COOPER_PAIR_HOPPING_IDS(Cooper_pair_hopping_names, Cooper_pair_hopping_energies, Cooper_pair_hoppings, n_cooper_pair_hoppings)
  IMPLICIT NONE
  INTEGER(INT32) :: n_cooper_pair_hoppings
  CHARACTER(LEN=100), INTENT(IN) :: Cooper_pair_hopping_names(n_cooper_pair_hoppings)
  REAL(real64), INTENT(IN) :: Cooper_pair_hopping_energies(n_cooper_pair_hoppings)
  TYPE(cooper_pair_hopping_t), INTENT(OUT) :: Cooper_pair_hoppings(n_cooper_pair_hoppings)

  INTEGER(INT32) :: i, cooper_pair_hopping_id_to_assign

  WRITE (log_string, *) "Assigning cooper pair hopping ids"
  LOG_INFO(log_string)

  DO i = 1, n_cooper_pair_hoppings
    cooper_pair_hopping_id_to_assign = 0
    SELECT CASE (TRIM(Cooper_pair_hopping_names(i)))
    CASE ("intraband")
      cooper_pair_hopping_id_to_assign = COOPER_PAIR_HOP_INTRABAND_ID
    CASE DEFAULT
      WRITE (log_string, *) "Unknown cooper pair hopping name: ", TRIM(Cooper_pair_hopping_names(i))
      LOG_ABNORMAL(log_string)
    END SELECT

    Cooper_pair_hoppings(i) = Cooper_pair_hopping_t(cooper_pair_hopping_id_to_assign, Cooper_pair_hopping_energies(i))
  END DO

END SUBROUTINE ASSIGN_COOPER_PAIR_HOPPING_IDS

SUBROUTINE ASSIGN_INTERACTION_IDS(Interaction_names, Interaction_energies, Interactions, n_interactions)
  IMPLICIT NONE
  INTEGER(INT32) :: n_interactions
  CHARACTER(LEN=100), INTENT(IN) :: Interaction_names(n_interactions)
  REAL(real64), INTENT(IN) :: Interaction_energies(n_interactions)
  TYPE(interaction_t), INTENT(OUT) :: Interactions(n_interactions)

  INTEGER(INT32) :: i, interaction_id_to_assign

  WRITE (log_string, *) "Assigning interaction ids"
  LOG_INFO(log_string)

  DO i = 1, n_interactions
    interaction_id_to_assign = 0
    SELECT CASE (TRIM(Interaction_names(i)))
    CASE ("intraband")
      interaction_id_to_assign = INT_INTRABAND_ID
    CASE ("interband")
      interaction_id_to_assign = INT_INTERBAND_ID
    CASE DEFAULT
      WRITE (log_string, *) "Unknown interaction name: ", TRIM(Interaction_names(i))
      LOG_ABNORMAL(log_string)
    END SELECT

    Interactions(i) = Interaction_t(interaction_id_to_assign, Interaction_energies(i))
  END DO

END SUBROUTINE ASSIGN_INTERACTION_IDS

SUBROUTINE CONSTRUCT_INTERACTION_TENSOR(Interactions, n_interactions, J_tensor, dim_tensor)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: n_interactions
  INTEGER(INT32), INTENT(IN) :: dim_tensor
  TYPE(interaction_t), INTENT(IN) :: Interactions(n_interactions)
  REAL(real64), INTENT(OUT) :: J_tensor(dim_tensor, dim_tensor, dim_tensor, dim_tensor)

  INTEGER(INT32) :: i, j
  WRITE (log_string, *) "Constructing interaction tensor"
  LOG_INFO(log_string)

  DO i = 1, n_interactions
    SELECT CASE (Interactions(i) % interaction_id)
    CASE (INT_INTRABAND_ID)
      CALL COMPUTE_INTRABAND_INTERACTION(Interactions(i) % energy, J_tensor, dim_tensor)
    CASE (INT_INTERBAND_ID)
      CALL COMPUTE_INTERBAND_INTERACTION(Interactions(i) % energy, J_tensor, dim_tensor)
    CASE DEFAULT
      WRITE (log_string, *) "Unknown interaction id: ", Interactions(i) % interaction_id
      LOG_ERROR(log_string)
    END SELECT
  END DO

  WRITE (log_string, *) "Computed all interactions"
  LOG_INFO(log_string)

END SUBROUTINE CONSTRUCT_INTERACTION_TENSOR

SUBROUTINE COMPUTE_COOPER_PAIR_HOPPINGS(Cooper_pair_hopping, n_cooper_pair_hoppings, Gamma_matrix, dim_gamma)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: n_cooper_pair_hoppings
  INTEGER(INT32), INTENT(IN) :: dim_gamma
  TYPE(cooper_pair_hopping_t), INTENT(IN) :: Cooper_pair_hopping(n_cooper_pair_hoppings)
  COMPLEX(REAL64), INTENT(INOUT) :: Gamma_matrix(dim_gamma, dim_gamma)

  COMPLEX(REAL64) :: Gamma_matrix_before_hoppings(dim_gamma, dim_gamma)
  INTEGER(INT32) :: n, i, j

  Gamma_matrix_before_hoppings = Gamma_matrix

  WRITE (log_string, *) "Computing cooper pair hoppings"
  LOG_INFO(log_string)

  DO n = 1, n_cooper_pair_hoppings
    SELECT CASE (Cooper_pair_hopping(n) % hopping_id)
    CASE (COOPER_PAIR_HOP_INTRABAND_ID)
      WRITE (log_string, *) "Computing intraband cooper pair hopping with multiplier = ", Cooper_pair_hopping(n) % multiplier
      LOG_INFO(log_string)
      DO i = 1, dim_gamma
        DO j = 1, dim_gamma
          IF (i == j) CYCLE
          Gamma_matrix(i, i) = Gamma_matrix(i, i) + Gamma_matrix_before_hoppings(j, j) * Cooper_pair_hopping(n) % multiplier
        END DO
      END DO
    CASE DEFAULT
      WRITE (log_string, *) "Unknown cooper pair hopping id: ", Cooper_pair_hopping(n) % hopping_id
      LOG_ERROR(log_string)
    END SELECT
  END DO

END SUBROUTINE COMPUTE_COOPER_PAIR_HOPPINGS

SUBROUTINE COMPUTE_INTRABAND_INTERACTION(energy, J_tensor, dim_tensor)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: dim_tensor
  REAL(real64), INTENT(IN) :: energy
  REAL(real64), INTENT(INOUT) :: J_tensor(dim_tensor, dim_tensor, dim_tensor, dim_tensor)

  INTEGER(INT32) :: i

  WRITE (log_string, *) "Computing intraband interaction with energy = ", energy
  LOG_INFO(log_string)

  DO i = 1, dim_tensor
    J_tensor(i, i, i, i) = J_tensor(i, i, i, i) + energy
  END DO

END SUBROUTINE COMPUTE_INTRABAND_INTERACTION

SUBROUTINE COMPUTE_INTERBAND_INTERACTION(energy, J_tensor, dim_tensor)
  IMPLICIT NONE
  INTEGER(INT32), INTENT(IN) :: dim_tensor
  REAL(real64), INTENT(IN) :: energy
  REAL(real64), INTENT(INOUT) :: J_tensor(dim_tensor, dim_tensor, dim_tensor, dim_tensor)

  INTEGER(INT32) :: i, j

  WRITE (log_string, *) "Computing interband interaction with energy = ", energy
  LOG_INFO(log_string)

  !TODO: Implement this type of interaction as needed
  ! DO i = 1, dim_tensor
  !   J_tensor(i, i, i, i) = J_tensor(i, i, i, i) + energy
  ! END DO

END SUBROUTINE COMPUTE_INTERBAND_INTERACTION

END MODULE interaction_factory
