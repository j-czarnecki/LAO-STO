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

MODULE parameters
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
IMPLICIT NONE
SAVE

!Complex numbers
COMPLEX(REAL64), PARAMETER :: imag = CMPLX(0., 1., KIND=REAL64)
COMPLEX(REAL64), PARAMETER :: c_zero = (0.0d0, 0.0d0) ! So that compiler does not complain about type mismatch between 0 and imag
COMPLEX(REAL64), PARAMETER :: c_one = (1.0d0, 0.0d0)

!Constants for hexagonal lattice and looping
INTEGER(INT32), PARAMETER :: N_NEIGHBOURS = 3
INTEGER(INT32), PARAMETER :: N_NEXT_NEIGHBOURS = 6
INTEGER(INT32), PARAMETER :: N_ALL_NEIGHBOURS = N_NEIGHBOURS + N_NEXT_NEIGHBOURS
INTEGER(INT32), PARAMETER :: SPINS = 2

!Unit conversion factors
REAL(REAL64), PARAMETER :: meV2au = 1./27211.
REAL(REAL64), PARAMETER :: nm2au = 1./0.05292
REAL(REAL64), PARAMETER :: T2au = 4.254382E-6

!Physical constants [a.u.]
REAL(REAL64), PARAMETER :: k_B = 8.617333262 * 1e-5 * 1e3 * meV2au
REAL(REAL64), PARAMETER :: PI = 4 * ATAN(1.0d0)
REAL(REAL64), PARAMETER :: muB = 0.5

!Constants for LAO/STO [a.u.]
REAL(REAL64), PARAMETER :: A_TILDE = SQRT(2./3.) * 0.3905 * nm2au !length

! Integration constants
REAL(REAL64), PARAMETER :: K1_MAX = 1. !Full Brillouin zone to integrate over
REAL(REAL64), PARAMETER :: K2_MAX = 1.
REAL(REAL64), PARAMETER :: KX_MAX = 4.0d0 * PI / (3.0d0 * SQRT(3.0d0)) !This defines maximum kx of hexagon corresponding to the first Brillouin Zone
REAL(REAL64), PARAMETER :: KY_MAX = 2.0d0 * PI / 3.0d0 !This defines maximum ky of hexagon corresponding to the first Brillouin Zone
REAL(REAL64), PARAMETER :: R_K_MAX = 4.0d0 * PI / (3.0d0 * SQRT(3.0d0)) !This defines the radius of circle that the first Brillouin Zone hexagon is inscribed in.
REAL(REAL64), PARAMETER :: JACOBIAN = 8 * PI**2 / (3.*SQRT(3.0d0))
INTEGER(INT32), PARAMETER :: N_BZ_SECTIONS = 6

!Orbital angular momentum matrices in [111] direction (z || [111])
COMPLEX(REAL64), PARAMETER :: L_x(3, 3) = TRANSPOSE(RESHAPE([c_zero, c_zero, imag, &
                                                             c_zero, c_zero, imag, &
                                                             -imag, -imag, c_zero], &
                                                            [3, 3])) / SQRT(2.0d0)
COMPLEX(REAL64), PARAMETER :: L_y(3, 3) = TRANSPOSE(RESHAPE([c_zero, -2 * imag, -imag, &
                                                             2 * imag, c_zero, imag, &
                                                             imag, -imag, c_zero], &
                                                            [3, 3])) / SQRT(6.0d0)
COMPLEX(REAL64), PARAMETER :: L_z(3, 3) = TRANSPOSE(RESHAPE([c_zero, -imag, imag, &
                                                             imag, c_zero, -imag, &
                                                             -imag, imag, c_zero], &
                                                            [3, 3])) / SQRT(3.0d0)

!Pauli matrices
COMPLEX(REAL64), PARAMETER :: Sigma_x(2, 2) = TRANSPOSE(RESHAPE([c_zero, c_one, &
                                                                 c_one, c_zero], &
                                                                [2, 2]))
COMPLEX(REAL64), PARAMETER :: Sigma_y(2, 2) = TRANSPOSE(RESHAPE([c_zero, -imag, &
                                                                 imag, c_zero], &
                                                                [2, 2]))
COMPLEX(REAL64), PARAMETER :: Sigma_z(2, 2) = TRANSPOSE(RESHAPE([c_one, c_zero, &
                                                                 c_zero, -c_one], &
                                                                [2, 2]))

END MODULE parameters
