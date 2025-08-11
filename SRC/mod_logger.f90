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

MODULE mod_logger
USE omp_lib
IMPLICIT NONE
SAVE
INTEGER*4, PRIVATE, PARAMETER :: LOGGER_UNIT = 66
INTEGER*4, PRIVATE, PARAMETER :: MAX_LOG_LEN = 2000
REAL*8, PRIVATE :: T_START, T_END
CHARACTER(LEN=MAX_LOG_LEN) :: log_string
!$omp threadprivate(log_string)

CHARACTER(LEN=*), parameter :: license_lines(*) = [ &
                               "Copyright (C) 2025 Julian Czarnecki", &
                               "", &
                               "This program is free software: you can redistribute it and/or modify", &
                               "it under the terms of the GNU General Public License as published by", &
                               "the Free Software Foundation, either version 3 of the License, or", &
                               "(at your option) any later version.", &
                               "", &
                               "This program is distributed in the hope that it will be useful,", &
                               "but WITHOUT ANY WARRANTY; without even the implied warranty of", &
                               "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", &
                               "GNU General Public License for more details.", &
                               "", &
                               "You should have received a copy of the GNU General Public License", &
                               "along with this program.  If not, see <https://www.gnu.org/licenses/>.", &
                               "", &
                               "If you use this code for scientific research, please cite:", &
                               "J. Czarnecki et. al.,", &
                               """Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface"",", &
                               "arXiv:2508.05075 (2025).", &
                               "https://arxiv.org/abs/2508.05075" &
                               ]

CONTAINS

RECURSIVE SUBROUTINE INIT_LOGGER(filename)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER*4 :: i
  IF (filename == "") THEN
    OPEN (unit=LOGGER_UNIT, FILE="./log.log", FORM="FORMATTED", ACTION="WRITE")
  ELSE
    OPEN (unit=LOGGER_UNIT, FILE="./"//TRIM(filename)//".log", FORM="FORMATTED", ACTION="WRITE")
  END IF
  DO i = 1, SIZE(license_lines)
    WRITE (LOGGER_UNIT, *) TRIM(license_lines(i))
  END DO
  WRITE (LOGGER_UNIT, *) "==== START ===="
  FLUSH (LOGGER_UNIT)
  CALL CPU_TIME(T_START)
END SUBROUTINE

RECURSIVE SUBROUTINE CLOSE_LOGGER()
  CALL CPU_TIME(T_END)
  WRITE (LOGGER_UNIT, *) "Time of simulation (seconds): ", T_END - T_START
  WRITE (LOGGER_UNIT, *) "==== END ===="
  CLOSE (LOGGER_UNIT)

END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_DEBUG(logMsg)
  CHARACTER(LEN=*), INTENT(IN) :: logMsg
  INTEGER*4 :: Values(8)
  CALL DATE_AND_TIME(VALUES=Values)
  !Add time printing
  WRITE (LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
  &'TID = ', omp_get_thread_num(), ' ',&
  &'DEBUG: '//TRIM(logMsg)
  FLUSH (LOGGER_UNIT)
END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_INFO(logMsg)
  CHARACTER(LEN=*), INTENT(IN) :: logMsg
  INTEGER*4 :: Values(8)

  CALL DATE_AND_TIME(VALUES=Values)
  !Add time printing
  WRITE (LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
  &'TID = ', omp_get_thread_num(), ' ',&
  &"INFO: "//TRIM(logMsg)
  FLUSH (LOGGER_UNIT)
END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_ABNORMAL(logMsg)
  CHARACTER(LEN=*), INTENT(IN) :: logMsg
  INTEGER*4 :: Values(8)

  CALL DATE_AND_TIME(VALUES=Values)
  !Add time printing
  WRITE (LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
  &'TID = ', omp_get_thread_num(), ' ',&
  &"ABNORMAL: "//TRIM(logMsg)
  FLUSH (LOGGER_UNIT)

END SUBROUTINE

RECURSIVE SUBROUTINE LOG_STRING_ERROR(logMsg)
  CHARACTER(LEN=*), INTENT(IN) :: logMsg
  INTEGER*4 :: Values(8)

  CALL DATE_AND_TIME(VALUES=Values)
  !Add time printing
  WRITE (LOGGER_UNIT, '(7(I0, a), a, I0, a, a)') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' ',&
  &'TID = ', omp_get_thread_num(), ' ',&
  &"ERROR: "//TRIM(logMsg)
  FLUSH (LOGGER_UNIT)

END SUBROUTINE

END MODULE
