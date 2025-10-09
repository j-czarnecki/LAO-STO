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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file contains macros definitions used in the code
! It should be included in every file that uses those macros
! via #include "macros_def.f90"
! It makes use of logger functions defined in logger.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Debug traces only if the DEBUG flag is active
! One has to specify -fpp flag during compilation to include preprocessor directives
! -DDEBUG to include debug traces
#ifdef DEBUG
#define LOG_DEBUG(logMsg) CALL LOG_STRING_DEBUG(logMsg)
#else
#define LOG_DEBUG(logMsg)
#endif

! Always define info traces
#define LOG_INFO(logMsg) CALL LOG_STRING_INFO(logMsg)
! Always define abnormal traces
#define LOG_ABNORMAL(logMsg) CALL LOG_STRING_ABNORMAL(logMsg)
! Always define error traces
#define LOG_ERROR(logMsg) CALL LOG_STRING_ERROR(logMsg)