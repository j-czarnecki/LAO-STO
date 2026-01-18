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
MODULE topology
use, intrinsic :: iso_fortran_env, only: real64, int8, int16, int32, int64
USE hamiltonians
USE parameters
USE utilities
USE writers
USE reader
USE local_integrand
USE self_consistency
USE logger
USE types
IMPLICIT NONE
CONTAINS

SUBROUTINE CALCULATE_CHERN_PARAMS(chern)
  !! Calculates Chern Params, based on https://arxiv.org/abs/cond-mat/0503172
  !! Not adapted to multiband systems
  TYPE(post_chern_number_t), INTENT(IN) :: chern
  ! INTEGER(INT32), INTENT(IN) :: Nk1 !! Number of divisions along k1
  ! INTEGER(INT32), INTENT(IN) :: Nk2 !! Number of divisions along k2
  !INTEGER(INT32), INTENT(IN) :: HamDim !! DImension of the hamiltonian to be diagonalized (e.g. 4 for simple hellical, 24 for LAO-STO)
  !CHARACTER(LEN=*) :: inputPath
  TYPE(sc_input_params_t) :: sc_input
  COMPLEX(REAL64), ALLOCATABLE :: Psi(:, :, :, :)
  COMPLEX(REAL64), ALLOCATABLE :: U1_chern(:, :), U2_chern(:, :), U3_chern(:, :), U4_chern(:, :)
  !COMPLEX(REAL64) :: U1_chern(DIM / 2, DIM / 2), U2_chern(DIM / 2, DIM / 2), U3_chern(DIM / 2, DIM / 2), U4_chern(DIM / 2, DIM / 2)
  INTEGER(INT32) :: i, j, a, b, m, n
  REAL(REAL64) :: potChem
  REAL(REAL64) :: Bfield(3)
  COMPLEX(REAL64) :: links

  COMPLEX(REAL64) :: f_12, det1, det2, det3, det4

  CALL GET_INPUT(TRIM(chern % path)//"input.nml", sc_input)

  ASSOCIATE (DIM_POSITIVE_K => sc_input % discretization % derived % DIM_POSITIVE_K, &
            & DIM => sc_input % discretization % derived % DIM)
    ALLOCATE (Psi(-chern % Nk_points / 2:chern % Nk_points / 2, 2, DIM, DIM))
    ALLOCATE (U1_chern(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (U2_chern(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (U3_chern(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE (U4_chern(DIM_POSITIVE_K, DIM_POSITIVE_K))
  END ASSOCIATE
  i = 0
  j = 0
  a = 0
  b = 0
  n = 0
  m = 0

  potChem = 0 * meV2au !Only for testing in hellical gap
  Bfield = (/0.0 * T2au, 0.0 * T2au, 0.0 * T2au/)

  !PRINT*, "Entered chern params"
  !PRINT*, "Ham dim ", DIM
  !PRINT*, "DIM/2", DIM/2
  !PRINT*, "Nk1/2", Nk1/2
  f_12 = CMPLX(0., 0., KIND=REAL64)
  !Calculate Chern numbers
  DO j = -chern % Nk_points / 2, chern % Nk_points / 2 - 1
    !This is for memory optimization. I dont have to keep all values of Psi over Brillouin Zone.
    !Instead I need values for current row and one row above:
    ! j = 0
    ! ********************************
    ! ********************************
    ! ********************************
    ! ********************************
    ! ********************************
    ! ********************************
    ! ******************************** <- this too
    ! ******************************** <- this I need
    !In next iteration I can forget about bottom row, the second one becomes the lower one
    ! and I have to calculate one above
    ! j = 1
    ! ********************************
    ! ********************************
    ! ********************************
    ! ********************************
    ! ********************************
    ! ******************************** <- this has to be calculated, Psi (:,2,:,:)
    ! ******************************** <- this becomes Psi(:,1,:,:)
    ! ******************************** <- this I can forget
    !It could be improved to keep Nk1 + 1 values, but for now I hope it is not necessary
    IF (j .EQ. 0) THEN
      DO n = -chern % Nk_points / 2, chern % Nk_points / 2
        CALL LAO_STO_CHERN_ENERGIES(chern % Nk_points, chern % Nk_points, n, j, chern % path, sc_input, Psi(n, 1, :, :)) !First row
        CALL LAO_STO_CHERN_ENERGIES(chern % Nk_points, chern % Nk_points, n, j + 1, chern % path, sc_input, Psi(n, 2, :, :)) !Second row

        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j, Psi(n,1,:,:))
        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j + 1, Psi(n,2,:,:))
      END DO
    ELSE
      Psi(:, 1, :, :) = Psi(:, 2, :, :)
      DO n = -chern % Nk_points / 2, chern % Nk_points / 2
        CALL LAO_STO_CHERN_ENERGIES(chern % Nk_points, chern % Nk_points, n, j + 1, chern % path, sc_input, Psi(n, 2, :, :)) !Next row
        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j+1, Psi(n,2,:,:)) !Next row
      END DO
    END IF

    DO i = -chern % Nk_points / 2, chern % Nk_points / 2 - 1
      !PRINT*, i, j
      !Calculate U matrices for chern numbers
      U1_chern = CMPLX(0.0, 0., KIND=REAL64)
      U2_chern = CMPLX(0.0, 0., KIND=REAL64)
      U3_chern = CMPLX(0.0, 0., KIND=REAL64)
      U4_chern = CMPLX(0.0, 0., KIND=REAL64)

      DO a = 1, sc_input % discretization % derived % DIM_POSITIVE_K
        DO b = 1, sc_input % discretization % derived % DIM_POSITIVE_K
          U1_chern(a, b) = SUM(CONJG(Psi(i, 1, :, a)) * Psi(i + 1, 1, :, b))
          U2_chern(a, b) = SUM(CONJG(Psi(i + 1, 1, :, a)) * Psi(i + 1, 2, :, b))
          U3_chern(a, b) = SUM(CONJG(Psi(i + 1, 2, :, a)) * Psi(i, 2, :, b))
          U4_chern(a, b) = SUM(CONJG(Psi(i, 2, :, a)) * Psi(i, 1, :, b))
        END DO
      END DO

      det1 = det(U1_chern(:, :), sc_input % discretization % derived % DIM_POSITIVE_K)
      IF (det1 .ne. 0.) THEN
        det1 = det1 / ABS(det1)
      END IF

      det2 = det(U2_chern(:, :), sc_input % discretization % derived % DIM_POSITIVE_K)
      IF (det2 .ne. 0.) THEN
        det2 = det2 / ABS(det2)
      END IF

      det3 = det(U3_chern(:, :), sc_input % discretization % derived % DIM_POSITIVE_K)
      IF (det3 .ne. 0.) THEN
        det3 = det3 / ABS(det3)
      END IF

      det4 = det(U4_chern(:, :), sc_input % discretization % derived % DIM_POSITIVE_K)
      IF (det4 .ne. 0.) THEN
        det4 = det4 / ABS(det4)
      END IF

      links = det1 * det2 * det3 * det4
      f_12 = f_12 + ATAN(AIMAG(links), REAL(links))

    END DO
  END DO

  OPEN (unit=9, FILE=TRIM(chern % path)//"OutputData/ChernNumber.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) f_12 / (2 * PI)
  CLOSE (9)
  !PRINT*, "Chern number is ", f_12/(2*PI)

  DEALLOCATE (Psi)
  DEALLOCATE (U1_chern)
  DEALLOCATE (U2_chern)
  DEALLOCATE (U3_chern)
  DEALLOCATE (U4_chern)
  IF (ALLOCATED(sc_input % physical % subband_params % V_layer)) DEALLOCATE (sc_input % physical % subband_params % V_layer)
  IF (ALLOCATED(sc_input % physical % subband_params % Subband_energies)) DEALLOCATE (sc_input % physical % subband_params % Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_CHERN_PARAMS

SUBROUTINE HELLICAL_TEST_CHERN(potChem, B, Nk1, Nk2, i, j, U_transformation)
  INTEGER(INT32), PARAMETER :: HamDim = 4
  REAL(REAL64), PARAMETER :: g = 5.0d0
  REAL(REAL64), PARAMETER :: tHop = 200 * meV2au
  REAL(REAL64), PARAMETER :: alphaSOC = 100 * meV2au
  REAL(REAL64), PARAMETER :: gammaSC = CMPLX(0.5 * meV2au, 0.0d0, KIND=REAL64)

  REAL(REAL64), INTENT(IN) :: potChem
  REAL(REAL64), INTENT(IN) :: B(3) ![Bx, By, Bz]
  INTEGER(INT32), INTENT(IN) :: Nk1, Nk2, i, j

  COMPLEX(REAL64), INTENT(INOUT) :: U_transformation(HamDim, HamDim)
  REAL(REAL64) :: Energies(HamDim)
  REAL(REAL64) :: dkx, dky
  REAL(REAL64) :: kx, ky
  INTEGER(INT32) :: m, n

  COMPLEX(REAL64) :: Hamiltonian(HamDim, HamDim)

  dkx = 2.0d0 * PI / Nk1
  dky = 2.0d0 * PI / Nk2

  Energies(:) = 0.0d0
  U_transformation(:, :) = CMPLX(0.0d0, 0.0d0, KIND=REAL64)
  Hamiltonian(:, :) = 0.0d0
  !Calculate all points in Brillouin zone
  kx = i * dkx
  ky = j * dky

  !Diagonal terms
  !Electrons H(k)
  Hamiltonian(1, 1) = 2.0 * tHop * (1.-COS(kx)) + 2.0 * tHop * (1 - COS(ky)) - potChem + 0.5 * muB * g * B(3)
  Hamiltonian(2, 2) = 2.0 * tHop * (1.-COS(kx)) + 2.0 * tHop * (1 - COS(ky)) - potChem - 0.5 * muB * g * B(3)

  !Holes -H*(-k)
  Hamiltonian(3, 3) = -(2.0 * tHop * (1.-COS(-kx)) + 2.0 * tHop * (1 - COS(-ky)) - potChem + 0.5 * muB * g * B(3))
  Hamiltonian(4, 4) = -(2.0 * tHop * (1.-COS(-kx)) + 2.0 * tHop * (1 - COS(-ky)) - potChem - 0.5 * muB * g * B(3))

  !Spin-orbit coupling
  Hamiltonian(1, 2) = 0.5 * mub * g * (B(1) - imag * B(2)) + alphaSOC * (SIN(kx) + imag * SIN(ky))
  Hamiltonian(3, 4) = -(0.5 * mub * g * (B(1) + imag * B(2)) + alphaSOC * (SIN(-kx) - imag * SIN(-ky)))

  !Superconductivity
  Hamiltonian(1, 4) = gammaSC
  Hamiltonian(2, 3) = -gammaSC

  Hamiltonian(:, :) = 0.5 * Hamiltonian(:, :)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, HamDim)

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:, :), Energies(:), U_transformation(:, :), HamDim)

  CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies(:), U_transformation(:, :), HamDim)

END SUBROUTINE HELLICAL_TEST_CHERN

SUBROUTINE LAO_STO_CHERN_ENERGIES(Nk1, Nk2, i, j, inputPath, sc_input, U_transformation)
    !! This subroutine calculates energies and wavefunctions of LAO-STO in [111] direction.
    !! Returns sorted wavefunctions in (i,j) point of the Brillouin zone.
  INTEGER(INT32), INTENT(IN) :: Nk1 !! Number of divisions of Brillouin zone in direction k1.
  INTEGER(INT32), INTENT(IN) :: Nk2 !! Number of divisions of Brillouin zone in direction k2.
  INTEGER(INT32), INTENT(IN) :: i !! Curent point k1_i in the Brillouin zone
  INTEGER(INT32), INTENT(IN) :: j !! Curent point k2_j in the Brillouin zone
  CHARACTER(LEN=*), INTENT(IN) :: inputPath !! Directory of the run, where input.nml should be placed.
                                            !! It contains material information and physical parameter of calculation:
                                            !! Fermi energy, temperature etc.

  TYPE(sc_input_params_t), INTENT(IN) :: sc_input
  COMPLEX(REAL64), INTENT(OUT) :: U_transformation(sc_input % discretization % derived % DIM, &
                                             & sc_input % discretization % derived % DIM) !! Matrix containing eigenvectors stored in consecutive columns.
                                                                                          !! On output sorted based on energies from lowest to highest.

  REAL(REAL64) :: Energies(sc_input % discretization % derived % DIM) !! Eigenvalues of the hamiltonian

  COMPLEX(REAL64) :: Hamiltonian(sc_input % discretization % derived % DIM,&
                                        & sc_input % discretization % derived % DIM)
  COMPLEX(REAL64) :: Hamiltonian_const(sc_input % discretization % derived % DIM,&
                                              & sc_input % discretization % derived % DIM)
  COMPLEX(REAL64) :: Gamma_SC(sc_input % discretization % ORBITALS,&
                                    & N_ALL_NEIGHBOURS,&
                                    & SPINS,&
                                    & SPINS,&
                                    & sc_input % discretization % derived % LAYER_COUPLINGS,&
                                    & sc_input % discretization % SUBBANDS)
  REAL(REAL64) :: Charge_dens(sc_input % discretization % derived % DIM_POSITIVE_K)
  REAL(REAL64) :: k1, k2, kx, ky
  REAL(REAL64) :: dk1_Chern, dk2_Chern
  INTEGER(INT32) :: n
  LOGICAL :: fileExists

  !PRINT*, "Allocation ended"
  dk1_Chern = K1_MAX / Nk1
  dk2_Chern = K2_MAX / Nk2

  U_transformation = CMPLX(0., 0., KIND=REAL64)
  Energies = 0.
  Gamma_SC = 0.
  Charge_dens = 0.

  !PRINT*, "Entered chern energies"
  !Get parameters from simulation

  CALL GET_SAFE_GAMMA_SC(Gamma_SC, inputPath, sc_input % discretization)
  CALL GET_SAFE_CHARGE_DENS(Charge_dens, inputPath, sc_input % discretization)

  !Gamma_SC(:,:,1,:) = 100.0d0 * meV2au
  !Gamma_SC(:,:,2,:) = -100.0d0 * meV2au

  !Computing k-independent terms
  Hamiltonian_const = CMPLX(0., 0., KIND=REAL64)
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const, sc_input % discretization, sc_input % physical)

  !Calculate eigenvalues and eigenvectors to later compute chern numbers
  k1 = i * dk1_Chern
  k2 = j * dk2_Chern

  kx = 2.*PI / (SQRT(3.0d0)) * k1
  ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
  Hamiltonian(:, :) = CMPLX(0., 0., KIND=REAL64)
  CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky, sc_input % discretization, sc_input % physical)
  CALL COMPUTE_HUBBARD(Hamiltonian, &
                      & Charge_dens, &
                      & sc_input % physical % subband_params % U_HUB, &
                      & sc_input % physical % subband_params % V_HUB, &
                      & sc_input % discretization)
  CALL COMPUTE_SC(Hamiltonian, kx, ky, Gamma_SC, sc_input % discretization)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, sc_input % discretization % derived % DIM) !This is not needed, since ZHEEV takes only upper triangle

  Hamiltonian = 0.5 * (Hamiltonian_const + Hamiltonian) !Should by multiplied by 0.5 if in Nambu space

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian, Energies, U_transformation, sc_input % discretization % derived % DIM)

  CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, U_transformation, sc_input % discretization % derived % DIM)

END SUBROUTINE LAO_STO_CHERN_ENERGIES

SUBROUTINE SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, Psi, HamDim)
  INTEGER(INT32), INTENT(IN) :: HamDim
  COMPLEX(REAL64), INTENT(INOUT) :: Psi(HamDim, HamDim)
  REAL(REAL64), INTENT(INOUT) :: Energies(HamDim)

  INTEGER(INT32) :: i, j
  REAL(REAL64) :: tmpEnergy
  COMPLEX(REAL64) :: tmpPsi(HamDim)

  DO i = 1, HamDim
    DO j = 1, HamDim - 1
    IF (Energies(j) .GT. Energies(j + 1)) THEN
      !Swap energies
      tmpEnergy = Energies(j)
      Energies(j) = Energies(j + 1)
      Energies(j + 1) = tmpEnergy

      !Swap wavefunctions
      tmpPsi(:) = Psi(:, j)
      Psi(:, j) = Psi(:, j + 1)
      Psi(:, j + 1) = tmpPsi
    END IF
    END DO
  END DO

  DO i = 1, HamDim
    Psi(:, i) = Psi(:, i) / SUM(ABS(Psi(:, i))**2)
  END DO

  ! DO i = 1, HamDim
  !     !PRINT*, "Psi(1,i) ", Psi(1,i)
  !     Psi(:,i) = Psi(:,i) / (Psi(1,i) / ABS(Psi(i,i)))
  ! END DO

END SUBROUTINE SORT_ENERGIES_AND_WAVEFUNCTIONS

END MODULE topology
