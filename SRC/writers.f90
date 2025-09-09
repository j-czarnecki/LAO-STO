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
USE parameters
USE reader
USe types
IMPLICIT NONE
CONTAINS

SUBROUTINE PRINT_HAMILTONIAN(Hamiltonian, N, name)
  INTEGER*4, INTENT(IN) :: N
  COMPLEX*16, INTENT(IN) :: Hamiltonian(N, N)
  CHARACTER(LEN=*), INTENT(IN) :: name
  CHARACTER(LEN=20) :: output_format
  INTEGER*4 :: i

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
  INTEGER*4, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: Energies(0:k1_steps, 0:k2_steps, N)
  REAL*8, INTENT(IN) :: dk1, dk2
  INTEGER*4, INTENT(IN) :: k1_steps, k2_steps
  REAL*8 :: k1, k2, kx, ky

  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format
  INTEGER*4 :: i, j, l

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
  COMPLEX*16, INTENT(IN) :: Gamma_SC(discretization % ORBITALS, &
                                    & N_ALL_NEIGHBOURS, &
                                    & SPINS, &
                                    & SPINS, &
                                    & discretization % derived % LAYER_COUPLINGS, &
                                    & discretization % SUBBANDS)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format

  INTEGER*4 :: orb, j, spin1, spin2, lat, band
  output_format = '(6I5, 2E15.5)'

  !Printing SC gammas in [meV]
  OPEN (unit=9, FILE="./OutputData/"//filename//".dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#band spin1 spin2 neighbour lattice orbital Re(Gamma) Im(Gamma)"
  DO band = 1, discretization % SUBBANDS
    DO spin1 = 1, SPINS
      DO spin2 = 1, SPINS
        !Do this, because we store LAYER_COUPLINGS for nearest neighbours, but only SUBLATTICES for next-to-nearest
        DO j = 1, N_NEIGHBOURS
          DO lat = 1, discretization % derived % LAYER_COUPLINGS
            DO orb = 1, discretization % ORBITALS
              WRITE (9, output_format) band, spin1, spin2, j, lat, orb, &
              & REAL(Gamma_SC(orb, j, spin1, spin2, lat, band)) / meV2au, &
              & AIMAG(Gamma_SC(orb, j, spin1, spin2, lat, band)) / meV2au
            END DO
          END DO
          WRITE (9, *)
          WRITE (9, *)
        END DO
        DO j = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
          DO lat = 1, discretization % SUBLATTICES
            DO orb = 1, discretization % ORBITALS
              WRITE (9, output_format) band, spin1, spin2, j, lat, orb, &
              & REAL(Gamma_SC(orb, j, spin1, spin2, lat, band)) / meV2au, &
              & AIMAG(Gamma_SC(orb, j, spin1, spin2, lat, band)) / meV2au
            END DO
          END DO
          WRITE (9, *)
          WRITE (9, *)
        END DO
      END DO
    END DO
  END DO
  CLOSE (9)
END SUBROUTINE PRINT_GAMMA

SUBROUTINE PRINT_CHARGE(Charge_dens, filename, discretization)
  TYPE(discretization_t), INTENT(IN) :: discretization
  REAL*8, INTENT(IN) :: Charge_dens(discretization % derived % DIM_POSITIVE_K, &
                                   & discretization % SUBBANDS)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=20) :: output_format
  INTEGER*4 :: spin, lat, orb, n, band

  output_format = '(4I5, 1E15.5)'
  OPEN (unit=9, FILE="./OutputData/"//filename//".dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) "#band spin lattice orbital Charge"
  DO band = 1, discretization % SUBBANDS
    n = 1
    DO spin = 1, 2
      DO lat = 1, discretization % SUBLATTICES
        DO orb = 1, discretization % ORBITALS
          WRITE (9, output_format) band, spin, lat, orb, Charge_dens(n, band)
          n = n + 1
        END DO
      END DO
    END DO
  END DO
  CLOSE (9)
END SUBROUTINE PRINT_CHARGE

END MODULE writers
