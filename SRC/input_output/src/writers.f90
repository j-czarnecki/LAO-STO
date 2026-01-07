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

MODULE writers
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE parameters
USE reader
USE types
USE utilities
IMPLICIT NONE
CONTAINS

SUBROUTINE PRINT_HAMILTONIAN(Hamiltonian, N, name)
  INTEGER(INT32), INTENT(IN) :: N
  COMPLEX(REAL64), INTENT(IN) :: Hamiltonian(N, N)
  CHARACTER(LEN=*), INTENT(IN) :: name
  CHARACTER(LEN=20) :: output_format
  INTEGER(INT32) :: i

  WRITE (output_format, '(A, I0, A)') '(', N, 'E15.5)'
  output_format = TRIM(output_format)

  OPEN (unit=9, FILE="./OutputData/"//TRIM(name)//"_real.dat", FORM="FORMATTED", ACTION="WRITE")
  OPEN (unit=10, FILE="./OutputData/"//TRIM(name)//"_imag.dat", FORM="FORMATTED", ACTION="WRITE")
  DO i = 1, N
    WRITE (9, output_format) REAL(Hamiltonian(i, :))
    WRITE (10, output_format) AIMAG(Hamiltonian(i, :))
  END DO
  CLOSE (9)
  CLOSE (10)

END SUBROUTINE PRINT_HAMILTONIAN

SUBROUTINE PRINT_ENERGIES(Energies, k1_steps, k2_steps, dk1, dk2, filename, N)
  INTEGER(INT32), INTENT(IN) :: N
  INTEGER(INT32), INTENT(IN) :: k1_steps, k2_steps
  REAL(REAL64), INTENT(IN) :: Energies(0:k1_steps, 0:k2_steps, N)
  REAL(REAL64), INTENT(IN) :: dk1, dk2
  REAL(REAL64) :: k1, k2, kx, ky

  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format
  INTEGER(INT32) :: i, j, l

  output_format = '(I5, 3E15.5)'

  OPEN (unit=9, FILE="./OutputData/"//filename//".dat", FORM="FORMATTED", ACTION="WRITE")
  DO l = 1, N
    DO i = 0, k1_steps
      DO j = 0, k2_steps
        k1 = i * dk1
        k2 = j * dk2

        kx = 2.*PI / (SQRT(3.0d0)) * k1
        ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
        WRITE (9, output_format) l, k1, k2, Energies(i, j, l) / meV2au
      END DO
      WRITE (9, *)
      WRITE (9, *)
    END DO
    WRITE (9, *)
    WRITE (9, *)
  END DO
  CLOSE (9)
END SUBROUTINE

SUBROUTINE PRINT_GAMMA(Gamma_SC, filename, discretization)
  TYPE(discretization_t), INTENT(IN) :: discretization
#ifndef BAND_BASIS
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(N_ALL_NEIGHBOURS + N_NEIGHBOURS, &
                                        & discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % SUBBANDS)
#else
  COMPLEX(REAL64), INTENT(IN) :: Gamma_SC(discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % derived % DIM_POSITIVE_K, &
                                        & discretization % SUBBANDS)
#endif
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format

  INTEGER(INT32) :: orb, j, spin1, spin2, lat, band
  INTEGER(INT32) :: i_band, j_band, neigh
  INTEGER(INT32) :: Degrees_of_freedom_row(4), Degrees_of_freedom_col(4)

  !Printing SC gammas in [meV]
  OPEN (unit=9, FILE="./OutputData/"//filename//".dat", FORM="FORMATTED", ACTION="WRITE")
#ifndef BAND_BASIS
  output_format = '(10I5, 2E15.5)'
  WRITE (9, '(100A)') "#band i_band j_band orb1 orb2 lat1 lat2 spin1 spin2 neighbour Re(Gamma) Im(Gamma)"
  DO band = 1, discretization % SUBBANDS
    Do i_band = 1, discretization % derived % DIM_POSITIVE_K
      Do j_band = 1, discretization % derived % DIM_POSITIVE_K
        Degrees_of_freedom_row = get_degress_of_freedom_from_index(i_band, discretization)
        Degrees_of_freedom_col = get_degress_of_freedom_from_index(j_band, discretization)
        DO neigh = 1, N_NEAREST_NEIGHBOURS + N_NEXT_NEIGHBOURS
          WRITE (9, output_format) band, &
          & i_band, j_band, &
          & Degrees_of_freedom_row(1), Degrees_of_freedom_col(1), &
          & Degrees_of_freedom_row(2), Degrees_of_freedom_col(2), &
          & Degrees_of_freedom_row(3), Degrees_of_freedom_col(3), &
          & neigh, &
          & REAL(Gamma_SC(neigh, i_band, j_band, band)) / meV2au, &
          & AIMAG(Gamma_SC(neigh, i_band, j_band, band)) / meV2au
        END DO
        WRITE (9, *)
        WRITE (9, *)
      END DO
    END DO
  END DO
#else
  output_format = '(3I5, 2E15.5)'
  WRITE (9, '(100A)') "#band i_band j_band Re(Gamma) Im(Gamma)"
  DO band = 1, discretization % SUBBANDS
    Do i_band = 1, discretization % derived % DIM_POSITIVE_K
      Do j_band = 1, discretization % derived % DIM_POSITIVE_K
        WRITE (9, output_format) band, &
        & i_band, j_band, &
        & REAL(Gamma_SC(i_band, j_band, band)) / meV2au, &
        & AIMAG(Gamma_SC(i_band, j_band, band)) / meV2au
      END DO
      WRITE (9, *)
      WRITE (9, *)
    END DO
  END DO
#endif
  CLOSE (9)
END SUBROUTINE PRINT_GAMMA

SUBROUTINE PRINT_CHARGE(Charge_dens, filename, discretization)
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL(REAL64), INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K, &
                                   & discretization % SUBBANDS)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format
  INTEGER(INT32) :: spin, lat, orb, n, band, i_band
  INTEGER(INT32) :: Degrees_of_freedom(4)

#ifndef BAND_BASIS
  output_format = '(5I5, 1E15.5)'
#else
  output_format = '(2I5, 1E15.5)'
#endif

  OPEN (unit=9, FILE="./OutputData/"//filename//".dat", FORM="FORMATTED", ACTION="WRITE")
#ifndef BAND_BASIS
  WRITE (9, '(100A)') "#band i_band orb lat spin Charge"
#else
  WRITE (9, '(100A)') "#band i_band Charge"
#endif

  DO band = 1, discretization % SUBBANDS
    DO i_band = 1, discretization % derived % DIM_POSITIVE_K
#ifndef BAND_BASIS
      Degrees_of_freedom = get_degress_of_freedom_from_index(i_band, discretization)
      WRITE (9, output_format) band, i_band, Degrees_of_freedom(1), Degrees_of_freedom(2), Degrees_of_freedom(3), Charge_dens(i_band, band)
#else
      WRITE (9, output_format) band, i_band, Charge_dens(i_band, band)
#endif
    END DO
  END DO
  CLOSE (9)
END SUBROUTINE PRINT_CHARGE

END MODULE writers
