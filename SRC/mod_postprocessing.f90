#include "macros_def.f90"
MODULE mod_postprocessing
USE mod_hamiltonians
USE mod_parameters
USE mod_utilities
USE mod_writers
USE mod_reader
USE mod_local_integrand
USE mod_self_consistency
USE mod_logger
IMPLICIT NONE
CONTAINS

SUBROUTINE CALCULATE_DOS(E_DOS_min, E_DOS_max, dE0, zeta_DOS, include_sc, Nk_points, N_refs, inputPath)

  !DOS calculation
  REAL*8, INTENT(IN) :: E_DOS_min, E_DOS_max, dE0
  REAL*8, INTENT(IN) :: zeta_DOS
  INTEGER*4, INTENT(IN) :: Nk_points, N_refs
  LOGICAL, INTENT(IN) :: include_sc
  CHARACTER(LEN=*), INTENT(IN) :: inputPath

  REAL*8 :: E0
  INTEGER*4 :: DOS_steps
  INTEGER*4 :: hamiltonian_dim

  COMPLEX*16, ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :), U_transformation(:, :)
  REAL*8, ALLOCATABLE :: Energies(:)

  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :)
  REAL*8, ALLOCATABLE :: Charge_dens(:, :)

  REAL*8, ALLOCATABLE :: DOS(:), DOS_local(:)
  CHARACTER(LEN=20) :: output_format

  REAL*8 :: k1, k2, kx, ky, dk1, dk2
  REAL*8 :: sc_multiplier
  INTEGER*4 :: i, j, k, n, lat, orb, orb_prime, spin, band

  INTEGER*4 :: i_ref, j_ref
  REAL*8 :: dk1_ref, dk2_ref

  INTEGER*4 :: points_within_energy_range, points_within_energy_range_local

  LOGICAL :: fileExists

  CALL GET_INPUT(TRIM(inputPath)//"input.nml")

  DOS_steps = INT((E_DOS_max - E_DOS_min) / dE0)
  dk1 = K1_MAX / Nk_points
  dk2 = K2_MAX / Nk_points

  dk1_ref = dk1 / (N_refs + 1)
  dk2_ref = dk2 / (N_refs + 1)

  !If superconductivity is to be included, we add Nambu space to the Hamiltonian and double the size.
  IF (include_sc) THEN
    hamiltonian_dim = DIM
    sc_multiplier = 0.5
  ELSE
    hamiltonian_dim = DIM_POSITIVE_K
    sc_multiplier = 1.0
  END IF

  ALLOCATE (Hamiltonian(DIM, DIM))
  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (Hamiltonian_const_band(DIM, DIM))
  ALLOCATE (U_transformation(hamiltonian_dim, hamiltonian_dim))
  ALLOCATE (Energies(hamiltonian_dim))
  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (DOS(0:DOS_steps))
  ALLOCATE (DOS_local(0:DOS_steps))

  Hamiltonian = DCMPLX(0., 0.)
  Hamiltonian_const = DCMPLX(0., 0.)
  Energies = 0.
  Gamma_SC = DCMPLX(0., 0.)
  Charge_dens = 0.
  DOS = 0.0d0
  DOS_local = 0.0d0
  INQUIRE (FILE=TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST=fileExists)
  IF (fileExists) THEN
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
  ELSE
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
  END IF

  IF (include_sc) THEN
    !Then we should also read Gamma_SC
    INQUIRE (FILE=TRIM(inputPath)//"OutputData/Gamma_SC_final.dat", EXIST=fileExists)
    IF (fileExists) THEN
      CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
    ELSE
      CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
    END IF
  END IF

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  DO band = 1, SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    !Adapt potential of given subband (energy difference due to quantization)
    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band)

    WRITE (log_string, *) "Calculating energies and integrating DOS..."
    LOG_INFO(log_string)
    !$omp parallel private(E0, k1, k2, kx, ky, Hamiltonian, Energies, U_transformation, DOS_local, points_within_energy_range_local)
    DOS_local = 0.0d0 !Initialize DOS_local for each thread!
    points_within_energy_range_local = 0
    !$omp do collapse(2)
    DO i = -Nk_points / 2, Nk_points / 2
      DO j = -Nk_points / 2, Nk_points / 2
        k1 = i * dk1
        k2 = j * dk2

        kx = 2.*PI / (SQRT(3.0d0)) * k1
        ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
        Hamiltonian(:, :) = DCMPLX(0., 0.)
        CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
        CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
        CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :, band))
        CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

        Hamiltonian(:, :) = sc_multiplier * (Hamiltonian_const_band + Hamiltonian) !Should by multiplied by 0.5 if in Nambu space

        CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), U_transformation(:, :), hamiltonian_dim)

        ! If lowest energy is beyond the range we are calculating the DOS for, skip
        ! This might be changed if magnetic field is to be introduced
        IF (MINVAL(ABS(Energies)) > E_DOS_max) CYCLE
        points_within_energy_range_local = points_within_energy_range_local + 1

        !Update DOS for current thread.
        DO n = 0, DOS_steps
          E0 = E_DOS_min + n * dE0
          DO k = 1, hamiltonian_dim
            DOS_local(n) = DOS_local(n) + dirac_delta(Energies(k), E0, zeta_DOS)
          END DO
        END DO

        !Add a grid refinement here, if lowest energy is in [E_DOS_min, E_DOS_max]
        !Center of the cell
        IF (i .lt. Nk_points / 2 .AND. j .lt. Nk_points / 2) THEN
          DO i_ref = 1, N_refs
            DO j_ref = 1, N_refs
              k1 = i * dk1 + i_ref * dk1_ref
              k2 = j * dk2 + j_ref * dk2_ref
              kx = 2.*PI / (SQRT(3.0d0)) * k1
              ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
              Hamiltonian(:, :) = DCMPLX(0., 0.)
              CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
              CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
              CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :, band))

              CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

              Hamiltonian(:, :) = sc_multiplier * (Hamiltonian_const_band + Hamiltonian) !Should by multiplied by 0.5 if in Nambu space

              CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), U_transformation(:, :), hamiltonian_dim)

              !Update DOS for current thread.
              DO n = 0, DOS_steps
                E0 = E_DOS_min + n * dE0
                DO k = 1, hamiltonian_dim
                  DOS_local(n) = DOS_local(n) + dirac_delta(Energies(k), E0, zeta_DOS)
                END DO
              END DO

            END DO
          END DO
        END IF

        IF ((i .lt. Nk_points / 2 .AND. j .lt. Nk_points / 2) .OR. (i .eq. Nk_points / 2)) THEN
          !Left edge without corner
          DO j_ref = 1, N_refs
            k1 = i * dk1
            k2 = j * dk2 + j_ref * dk2_ref
            kx = 2.*PI / (SQRT(3.0d0)) * k1
            ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
            Hamiltonian(:, :) = DCMPLX(0., 0.)
            CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
            CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
            CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :, band))

            CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

            Hamiltonian(:, :) = sc_multiplier * (Hamiltonian_const_band + Hamiltonian) !Should by multiplied by 0.5 if in Nambu space

            CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), U_transformation(:, :), hamiltonian_dim)

            !Update DOS for current thread.
            DO n = 0, DOS_steps
              E0 = E_DOS_min + n * dE0
              DO k = 1, hamiltonian_dim
                DOS_local(n) = DOS_local(n) + dirac_delta(Energies(k), E0, zeta_DOS)
              END DO
            END DO

          END DO
        END IF

        IF ((i .lt. Nk_points / 2 .AND. j .lt. Nk_points / 2) .OR. (j .eq. Nk_points / 2)) THEN
          !Bottom edge without corner
          DO i_ref = 1, N_refs
            k1 = i * dk1 + i_ref * dk1_ref
            k2 = j * dk2
            kx = 2.*PI / (SQRT(3.0d0)) * k1
            ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
            Hamiltonian(:, :) = DCMPLX(0., 0.)
            CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
            CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
            CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :, band))

            CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

            Hamiltonian(:, :) = sc_multiplier * (Hamiltonian_const_band + Hamiltonian) !Should by multiplied by 0.5 if in Nambu space

            CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(:), U_transformation(:, :), hamiltonian_dim)

            !Update DOS for current thread.
            DO n = 0, DOS_steps
              E0 = E_DOS_min + n * dE0
              DO k = 1, hamiltonian_dim
                DOS_local(n) = DOS_local(n) + dirac_delta(Energies(k), E0, zeta_DOS)
              END DO
            END DO

          END DO
        END IF

      END DO
    END DO
    !$omp end do

    !$omp critical (accumulate_dos)
    DO n = 0, DOS_steps
      DOS(n) = DOS(n) + DOS_local(n)
    END DO
    points_within_energy_range = points_within_energy_range + points_within_energy_range_local
    !$omp end critical (accumulate_dos)

    !$omp end parallel
  END DO

  ! Normalize the DOS
  IF (MAXVAL(DOS) .NE. 0.0d0) DOS = DOS / MAXVAL(DOS)

  WRITE (log_string, *) "Included points: ", points_within_energy_range, "out of ", Nk_points**2
  LOG_INFO(log_string)

  WRITE (log_string, *) "Writing DOS to file"
  LOG_INFO(log_string)

  output_format = '(2E15.5)'
  OPEN (unit=9, FILE=TRIM(inputPath)//"OutputData/DOS.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, '(A)') "#E[meV] DOS[a.u]"
  DO n = 0, DOS_steps
    E0 = E_DOS_min + n * dE0
    WRITE (9, output_format) E0 / meV2au, DOS(n)
  END DO
  CLOSE (9)

  DEALLOCATE (Hamiltonian)
  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (Hamiltonian_const_band)
  DEALLOCATE (U_transformation)
  DEALLOCATE (Energies)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (DOS)
  DEALLOCATE (DOS_local)
  IF (ALLOCATED(V_layer)) DEALLOCATE (V_layer)
  IF (ALLOCATED(Subband_energies)) DEALLOCATE (Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_DOS

SUBROUTINE CALCULATE_DISPERSION(inputPath, Nk_points, include_sc)
    !! Calculates dispersion relation in the first Brillouin zone.
    !! Takes physical parameters from input.nml from a directory specified by inputPath.
  CHARACTER(LEN=*), INTENT(IN) :: inputPath
  INTEGER*4, INTENT(IN) :: Nk_points
  LOGICAL, INTENT(IN) :: include_sc
  CHARACTER(LEN=20) :: output_format

  COMPLEX*16, ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :), U_transformation(:, :)
  REAL*8, ALLOCATABLE :: Energies(:, :, :)
  REAL*8, ALLOCATABLE :: Probability(:, :, :, :) ! |Psi^2|

  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :)
  REAL*8, ALLOCATABLE :: Charge_dens(:, :)

  REAL*8 :: k1, k2, kx, ky, dkx, dky
  REAL*8 :: sc_multiplier
  INTEGER*4 :: i, j, k, n, lat, orb, orb_prime, spin, l, m, band
  INTEGER*4 :: kx_steps, ky_steps
  REAL*8 :: yz_contribution, zx_contribution, xy_contribution
  REAL*8 :: lat1_contribution, lat2_contribution
  REAL*8, ALLOCATABLE :: Lat_contributions(:)
  REAL*8 :: spin_up_contribution, spin_down_contribution
  REAL*8 :: electron_contribution, hole_contribution
  REAL*8 :: brillouinZoneVertices(6, 2)

  LOGICAL :: fileExists
  INTEGER*4 :: hamiltonian_dim

  brillouinZoneVertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  brillouinZoneVertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  CALL GET_INPUT(TRIM(inputPath)//"input.nml")

  !If superconductivity is to be included, we add Nambu space to the Hamiltonian and double the size.
  IF (include_sc) THEN
    hamiltonian_dim = DIM
    sc_multiplier = 0.5d0
  ELSE
    hamiltonian_dim = DIM_POSITIVE_K
    sc_multiplier = 1.0d0
  END IF

  yz_contribution = 0.
  zx_contribution = 0.
  xy_contribution = 0.
  lat1_contribution = 0.
  lat2_contribution = 0.
  spin_up_contribution = 0.
  spin_down_contribution = 0.
  electron_contribution = 0.
  hole_contribution = 0.

  dkx = KX_MAX / Nk_points
  dky = KY_MAX / Nk_points

  kx_steps = INT(Nk_points)
  ky_steps = INT(Nk_points)

  WRITE (output_format, '(A, I0, A)') '(I5, ', 10 + SUBLATTICES, 'E15.5)'

  ALLOCATE (Hamiltonian(DIM, DIM))
  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (Hamiltonian_const_band(DIM, DIM))
  ALLOCATE (U_transformation(hamiltonian_dim, hamiltonian_dim))
  ALLOCATE (Probability(-kx_steps:kx_steps, -ky_steps:ky_steps, hamiltonian_dim, hamiltonian_dim))
  ALLOCATE (Energies(-kx_steps:kx_steps, -ky_steps:ky_steps, hamiltonian_dim))
  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Lat_contributions(SUBLATTICES))

  Hamiltonian(:, :) = DCMPLX(0., 0.)
  Hamiltonian_const(:, :) = DCMPLX(0., 0.)
  Probability(:, :, :, :) = 0.
  Energies(:, :, :) = 0.
  Gamma_SC = DCMPLX(0., 0.) * meV2au
  Charge_dens = 0.

  IF ((U_HUB .NE. 0.0) .OR. (V_HUB .NE. 0.0)) THEN
    INQUIRE (FILE=TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST=fileExists)
    IF (fileExists) THEN
      CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
    ELSE
      CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
    END IF
  END IF

  IF (include_sc) THEN
    !Then we should also read Gamma_SC
    INQUIRE (FILE=TRIM(inputPath)//"OutputData/Gamma_SC__final.dat", EXIST=fileExists)
    IF (fileExists) THEN
      CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
    ELSE
      CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
    END IF
  END IF

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  OPEN (unit=9, FILE=TRIM(inputPath)//"OutputData/Energies.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, '(A)') "#N kx[1/a] ky[1/a] Energy[meV] P(yz) P(zx) P(xy) P(lat1) P(lat2) ... P(latN) P(s_up) P(s_down) P(electron) P(hole)"

  DO band = 1, SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    !Adapt potential of given subband (energy difference due to quantization)
    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band)

    !$omp parallel private(kx, ky, Hamiltonian)
    !$omp do
    DO i = -kx_steps, kx_steps
      DO j = -ky_steps, ky_steps
        kx = i * dkx !* (2. * PI * 2./3.)
        ky = j * dky !* (2. * PI * 2./3.)
        IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

          Hamiltonian(:, :) = DCMPLX(0., 0.)
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky)
          CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens)
          CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC)
          CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle
          Hamiltonian(:, :) = sc_multiplier * (Hamiltonian_const_band(:, :) + Hamiltonian(:, :)) !Should by multiplied by 0.5 if in Nambu space

          !CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:DIM_POSITIVE_K,:DIM_POSITIVE_K), Energies(i,j,:), U_transformation(:,:), DIM_POSITIVE_K)
          !Probability(i,j,:,:) = ABS(U_transformation)**2

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim), Energies(i, j, :), hamiltonian_dim)
          Probability(i, j, :, :) = ABS(Hamiltonian(:hamiltonian_dim, :hamiltonian_dim))**2
        END IF
      END DO
    END DO
    !$omp end do
    !$omp end parallel

    DO l = 1, hamiltonian_dim
      DO i = -kx_steps, kx_steps
        DO j = -ky_steps, ky_steps
          kx = i * dkx !* (2. * PI * 2./3.)
          ky = j * dky !* (2. * PI * 2./3.)
          IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

            !Calculate specific contributions
            !Distinguishing orbital contributions
            yz_contribution = 0.
            zx_contribution = 0.
            xy_contribution = 0.
            DO n = 1, hamiltonian_dim, ORBITALS
              yz_contribution = yz_contribution + Probability(i, j, n, l)
              zx_contribution = zx_contribution + Probability(i, j, n + 1, l)
              xy_contribution = xy_contribution + Probability(i, j, n + 2, l)
            END DO

            !Distinguishing lattice contributions
            lat1_contribution = 0.
            lat2_contribution = 0.
            Lat_contributions(:) = 0.
            DO m = 0, SUBLATTICES - 1
              DO spin = 0, 1
                DO n = 1, ORBITALS
                  Lat_contributions(m + 1) = Lat_contributions(m + 1) + Probability(i, j, spin * TBA_DIM + m * ORBITALS + n, l)
                  IF (include_sc) THEN
                    Lat_contributions(m + 1) = Lat_contributions(m + 1) + Probability(i, j, DIM_POSITIVE_K + spin * TBA_DIM + m * ORBITALS + n, l)
                  END IF
                END DO
              END DO
            END DO

            !Distinguishing spin contributions
            spin_up_contribution = 0.
            spin_down_contribution = 0.
            DO n = 1, TBA_DIM
              spin_up_contribution = spin_up_contribution + Probability(i, j, n, l)
              spin_down_contribution = spin_down_contribution + Probability(i, j, TBA_DIM + n, l)
              IF (include_sc) THEN
                spin_up_contribution = spin_up_contribution + Probability(i, j, DIM_POSITIVE_K + n, l)
                spin_down_contribution = spin_down_contribution + Probability(i, j, DIM_POSITIVE_K + TBA_DIM + n, l)
              END IF
            END DO

            electron_contribution = 0.
            hole_contribution = 0.
            IF (include_sc) THEN
              DO n = 1, DIM_POSITIVE_K
                electron_contribution = electron_contribution + Probability(i, j, n, l)
                hole_contribution = hole_contribution + Probability(i, j, DIM_POSITIVE_K + n, l)
              END DO
            ELSE
              electron_contribution = 1.
              hole_contribution = 0.
            END IF

            WRITE (9, output_format) band * hamiltonian_dim + l, kx, ky, Energies(i, j, l) / meV2au, &
            & yz_contribution, zx_contribution, xy_contribution, &
            & (Lat_contributions(lat), lat=1, SUBLATTICES), &
            & spin_up_contribution, spin_down_contribution, &
            & electron_contribution, hole_contribution
          END IF
        END DO
      END DO
      ! WRITE(9,*)
      ! WRITE(9,*)
    END DO
  END DO !End of iteration over subbands
  CLOSE (9)

  DEALLOCATE (Hamiltonian)
  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (U_transformation)
  DEALLOCATE (Probability)
  DEALLOCATE (Energies)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (Lat_contributions)
  IF (ALLOCATED(V_layer)) DEALLOCATE (V_layer)
  IF (ALLOCATED(Subband_energies)) DEALLOCATE (Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_DISPERSION

SUBROUTINE CALCULATE_CHERN_PARAMS(Nk1, Nk2, inputPath)
    !! Calculates Chern Params, based on https://arxiv.org/abs/cond-mat/0503172
  INTEGER*4, INTENT(IN) :: Nk1 !! Number of divisions along k1
  INTEGER*4, INTENT(IN) :: Nk2 !! Number of divisions along k2
  !INTEGER*4, INTENT(IN) :: HamDim !! DImension of the hamiltonian to be diagonalized (e.g. 4 for simple hellical, 24 for LAO-STO)
  CHARACTER(LEN=*) :: inputPath

  COMPLEX*16, ALLOCATABLE :: Psi(:, :, :, :)
  COMPLEX*16 :: U1_chern(DIM / 2, DIM / 2), U2_chern(DIM / 2, DIM / 2), U3_chern(DIM / 2, DIM / 2), U4_chern(DIM / 2, DIM / 2)
  INTEGER*4 :: i, j, a, b, m, n
  REAL*8 :: potChem
  REAL*8 :: Bfield(3)
  COMPLEX*16 :: links

  COMPLEX*16 :: f_12, det1, det2, det3, det4

  ALLOCATE (Psi(-Nk1 / 2:Nk1 / 2, 2, DIM, DIM))
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
  f_12 = DCMPLX(0., 0.)
  !Calculate Chern numbers
  DO j = -Nk2 / 2, Nk2 / 2 - 1
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
      DO n = -Nk1 / 2, Nk1 / 2
        CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2, n, j, inputPath, Psi(n, 1, :, :)) !First row
        CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2, n, j + 1, inputPath, Psi(n, 2, :, :)) !Second row

        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j, Psi(n,1,:,:))
        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j + 1, Psi(n,2,:,:))
      END DO
    ELSE
      Psi(:, 1, :, :) = Psi(:, 2, :, :)
      DO n = -Nk1 / 2, Nk1 / 2
        CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2, n, j + 1, inputPath, Psi(n, 2, :, :)) !Next row
        !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j+1, Psi(n,2,:,:)) !Next row
      END DO
    END IF

    DO i = -Nk1 / 2, Nk1 / 2 - 1
      !PRINT*, i, j
      !Calculate U matrices for chern numbers
      U1_chern = DCMPLX(0.0, 0.)
      U2_chern = DCMPLX(0.0, 0.)
      U3_chern = DCMPLX(0.0, 0.)
      U4_chern = DCMPLX(0.0, 0.)

      DO a = 1, DIM / 2
        DO b = 1, DIM / 2
          U1_chern(a, b) = SUM(CONJG(Psi(i, 1, :, a)) * Psi(i + 1, 1, :, b))
          U2_chern(a, b) = SUM(CONJG(Psi(i + 1, 1, :, a)) * Psi(i + 1, 2, :, b))
          U3_chern(a, b) = SUM(CONJG(Psi(i + 1, 2, :, a)) * Psi(i, 2, :, b))
          U4_chern(a, b) = SUM(CONJG(Psi(i, 2, :, a)) * Psi(i, 1, :, b))
        END DO
      END DO

      det1 = det(U1_chern(:, :), DIM / 2)
      IF (det1 .ne. 0.) THEN
        det1 = det1 / ABS(det1)
      END IF

      det2 = det(U2_chern(:, :), DIM / 2)
      IF (det2 .ne. 0.) THEN
        det2 = det2 / ABS(det2)
      END IF

      det3 = det(U3_chern(:, :), DIM / 2)
      IF (det3 .ne. 0.) THEN
        det3 = det3 / ABS(det3)
      END IF

      det4 = det(U4_chern(:, :), DIM / 2)
      IF (det4 .ne. 0.) THEN
        det4 = det4 / ABS(det4)
      END IF

      links = det1 * det2 * det3 * det4
      f_12 = f_12 + ATAN(AIMAG(links), REAL(links))

    END DO
  END DO

  OPEN (unit=9, FILE=TRIM(inputPath)//"OutputData/ChernNumber.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) f_12 / (2 * PI)
  CLOSE (9)
  !PRINT*, "Chern number is ", f_12/(2*PI)

  DEALLOCATE (Psi)

END SUBROUTINE CALCULATE_CHERN_PARAMS

SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP(inputPath, dE, nBrillouinPoints)

  CHARACTER(LEN=*), INTENT(IN) :: inputPath !! This should be a path to folder where input.nml resides
  REAL*8, INTENT(IN) :: dE
  INTEGER*4, INTENT(IN) :: nBrillouinPoints
  CHARACTER(LEN=20) :: output_format

  COMPLEX*16, ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Hamiltonian_const_band(:, :)
  REAL*8, ALLOCATABLE :: Energies(:)

  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :)
  REAL*8, ALLOCATABLE :: Charge_dens(:, :)

  INTEGER*1, ALLOCATABLE :: IsFermiSurface(:, :) !! This indicates whether given (kx,ky) point is at Fermi surface
  INTEGER*2, ALLOCATABLE :: OrbitalAtFermiSurface(:, :) !! This indicates which state consitutes to the Fermi surface.
                                                         !! Because ZHEEV is used, energies are sorted from lowest to highest.
                                                         !! This implies that we do not recognize which spin, orbital etc.
                                                         !! constitutes to the Fermi surface, only "lowest", "second lowest" and so on.

  REAL*8 :: brillouinZoneVertices(6, 2)

  REAL*8 :: k1, k2, kx, ky, dkx, dky
  INTEGER*4 :: i, j, n, m, band
  INTEGER*4 :: kx_steps, ky_steps

  LOGICAL :: fileExists

  brillouinZoneVertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  brillouinZoneVertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  dkx = KX_MAX / nBrillouinPoints
  dky = KY_MAX / nBrillouinPoints

  kx_steps = INT(nBrillouinPoints)
  ky_steps = INT(nBrillouinPoints)

  !Get parameters from simulation
  CALL GET_INPUT(TRIM(inputPath)//"input.nml")

  ALLOCATE (Hamiltonian(DIM, DIM))
  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (Hamiltonian_const_band(DIM, DIM))
  ALLOCATE (Energies(DIM))
  ALLOCATE (IsFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
  ALLOCATE (OrbitalAtFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))

  Hamiltonian(:, :) = DCMPLX(0., 0.)
  Hamiltonian_const(:, :) = DCMPLX(0., 0.)
  Energies(:) = 0.
  Gamma_SC = DCMPLX(0., 0.) * meV2au
  Charge_dens = 0.

  output_format = '(3E15.5, I10)'

  !Calculation of superconducting gap
  INQUIRE (FILE=TRIM(inputPath)//"OutputData/Gamma_SC_final.dat", EXIST=fileExists)
  IF (fileExists) THEN
    CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
  ELSE
    CALL GET_GAMMA_SC(Gamma_SC, TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
  END IF

  INQUIRE (FILE=TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST=fileExists)
  IF (fileExists) THEN
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
  ELSE
    CALL GET_CHARGE_DENS(Charge_dens, TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
  END IF
  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  OPEN (unit=9, FILE=TRIM(inputPath)//"OutputData/SuperconductingGap.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (9, *) '#kx[1/a] ky[1/a] gap_SC[meV] N_orbital'
  DO band = 1, SUBBANDS
    WRITE (log_string, *) "Band: ", band
    LOG_INFO(log_string)

    IsFermiSurface(:, :) = 0
    OrbitalAtFermiSurface(:, :) = 0

    Hamiltonian_const_band = Hamiltonian_const
    CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band)

    !Dispersion relation in a normal state
    !$omp parallel private(kx, ky, Hamiltonian, Energies)
    !$omp do
    DO i = -kx_steps, kx_steps
      DO j = -ky_steps, ky_steps
        kx = i * dkx !* (2. * PI * 2./3.)
        ky = j * dky !* (2. * PI * 2./3.)

        IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

          Hamiltonian(:, :) = DCMPLX(0., 0.)
          Energies(:) = 0.
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
          CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
          !CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))

          CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

          Hamiltonian(:, :) = Hamiltonian_const_band(:, :) + Hamiltonian(:, :) !Should by multiplied by 0.5 if in Nambu space

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:DIM_POSITIVE_K, :DIM_POSITIVE_K), Energies(:DIM_POSITIVE_K), DIM_POSITIVE_K)
          !CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:DIM_POSITIVE_K, :DIM_POSITIVE_K), Energies(:DIM_POSITIVE_K), U_transformation(:DIM_POSITIVE_K, :DIM_POSITIVE_K), DIM_POSITIVE_K)

          !Check whether current wavevector is in the Fermi surface
          IF (MINVAL(ABS(Energies(:DIM_POSITIVE_K))) < dE) THEN
            IsFermiSurface(i, j) = 1
            OrbitalAtFermiSurface(i, j) = (band - 1) * DIM_POSITIVE_K + MINLOC(ABS(Energies(:DIM_POSITIVE_K)), 1)
          END IF
        END IF
      END DO
    END DO
    !$omp end do
    !$omp end parallel

    WRITE (log_string, *) "Fermi surface done"
    LOG_INFO(log_string)

    !$omp parallel private(kx, ky, Hamiltonian, Energies)
    !$omp do
    DO i = -kx_steps, kx_steps
      DO j = -ky_steps, ky_steps
        IF (IsFermiSurface(i, j) == 1) THEN
          kx = i * dkx
          ky = j * dky

          Hamiltonian(:, :) = DCMPLX(0., 0.)
          Energies(:) = 0.
          CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian(:, :), kx, ky)
          CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:, band))
          CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :, band))
          ! DO n = 1, DIM_POSITIVE_K
          !     IF (n .le. DIM_POSITIVE_K - ORBITALS) THEN
          !         Hamiltonian(n, ORBITALS + DIM_POSITIVE_K + n) = 20 * meV2au
          !     END IF
          !     IF (n .ge. ORBITALS) THEN
          !         Hamiltonian(n, DIM_POSITIVE_K + n - ORBITALS) = 20 * meV2au
          !     END IF
          ! END DO
          CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

          Hamiltonian(:, :) = 0.5 * (Hamiltonian_const_band(:, :) + Hamiltonian(:, :)) !Should by multiplied by 0.5 if in Nambu space

          CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:, :), Energies(:), DIM)
          !Write superconducting gap
          WRITE (9, output_format) kx, ky, ABS(Energies(DIM_POSITIVE_K) - Energies(DIM_POSITIVE_K + 1)) / meV2au, OrbitalAtFermiSurface(i, j)
        END IF
      END DO
    END DO
    !$omp end do
    !$omp end parallel
  END DO
  CLOSE (9)

  IF (ALLOCATED(V_layer)) DEALLOCATE (V_layer)
  IF (ALLOCATED(Subband_energies)) DEALLOCATE (Subband_energies) !Deallocate global variable

END SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP

SUBROUTINE CALCULATE_GAMMA_K(input_path, n_brillouin_points)
  CHARACTER(LEN=*), INTENT(IN) :: input_path !! This should be a path to folder where input.nml resides
  INTEGER*4, INTENT(IN) :: n_brillouin_points !! Number of steps taken in k-space.
                                                !! Integration is over interval -KX_MAX < ----- n_brillouin_points ----- 0 ----- n_brillouin_points ----- > KX_MAX
                                                !! So effectively 2N + 1 steps are taken in each direction

  COMPLEX*16, ALLOCATABLE :: Hamiltonian_const(:, :) !! k-indepndent and band-independent part of the Hamiltonian
  COMPLEX*16, ALLOCATABLE :: Hamiltonian_const_band(:, :) !! k-independent, band-dependent Hamiltonian
  COMPLEX*16, ALLOCATABLE :: Hamiltonian_dummy(:, :) !! This is only used to get matrix elements gamma that show up in the Hamiltonian
  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :) !! Supercondicting pairings read from a simulation
  COMPLEX*16, ALLOCATABLE :: Gamma_K(:, :, :, :, :) !! Superconducting pairing, determined at given k point
  REAL*8, ALLOCATABLE :: Charge_dens(:, :) !! Charge density read from a simulation
  COMPLEX*16, ALLOCATABLE :: Delta_local(:, :, :, :, :) !! Delta (pairing amplitudes) for given k point, integrand
  REAL*8, ALLOCATABLE :: Charge_dens_local(:, :) !! Charge density for given k point, integrand

  REAL*8 :: k1, k2, kx, ky, dkx, dky
  INTEGER*4 :: i, j, n, m, band
  INTEGER*4 :: kx_steps, ky_steps
  INTEGER*4 :: orb, neigh, spin, layer, lat, file_count
  INTEGER*4 :: orb_prime, band_prime
  INTEGER*4 :: row, col, row_inverse, col_inverse, row_nnn, col_nnn
  INTEGER*4 :: gamma_lat_index, gamma_spin_index

  COMPLEX*16 :: gamma_nn_12, gamma_nn_21, gamma_nnn

  CHARACTER(LEN=200) :: filename
  COMPLEX*16, ALLOCATABLE :: File_unit_mapping(:, :, :, :)

  LOGICAL :: file_exists

  REAL*8 :: Brillouin_zone_vertices(6, 2)

  Brillouin_zone_vertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  Brillouin_zone_vertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  dkx = 1 * KX_MAX / n_brillouin_points
  dky = 1 * KY_MAX / n_brillouin_points

  kx_steps = INT(n_brillouin_points)
  ky_steps = INT(n_brillouin_points)

  CALL GET_INPUT(TRIM(input_path)//"input.nml")

  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (Hamiltonian_const_band(DIM, DIM))
  ALLOCATE (Hamiltonian_dummy(DIM, DIM))
  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Gamma_K(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (Delta_local(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens_local(DIM_POSITIVE_K, SUBBANDS))
  ALLOCATE (File_unit_mapping(ORBITALS, 2, LAYER_COUPLINGS, SUBBANDS))

  Hamiltonian_const = DCMPLX(0., 0.)
  Hamiltonian_const_band = DCMPLX(0., 0.)
  Gamma_SC = DCMPLX(0., 0.)
  Charge_dens = 0.

  CALL GET_SAFE_GAMMA_SC(Gamma_SC, input_path)
  CALL GET_SAFE_CHARGE_DENS(Charge_dens, input_path)

  !Computing k-independent terms
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  !Opening all files I will write gammas to and create a mapping of file units to appropriate names
  file_count = 10
  DO orb = 1, ORBITALS
    DO spin = 1, 2
      DO layer = 1, LAYER_COUPLINGS
        DO band = 1, SUBBANDS
          File_unit_mapping(orb, spin, layer, band) = file_count
          WRITE (filename, "(4(A, I0))") "GammaK_orb", orb, "_spin", spin, "_layer", layer, "_band", band
          OPEN (unit=file_count, FILE=TRIM(input_path)//"OutputData/"//TRIM(filename)//".dat", FORM="FORMATTED", ACTION="WRITE")
          WRITE (file_count, '(A)') "#kx[1/a]   ky[1/a]     Re(Gamma_ham_nn)[meV]     Im(Gamma_ham_nn)[meV]    OPTIONAL(if layer <= SUBLATTICES)[Re(Gamma_ham_nnn)[meV]     Im(Gamma_ham_nnn)[meV]]     Re(Gamma_neigh1)[meV]   Im(Gamma_neigh1)[meV]   Re(Gamma_neigh2)[meV] ..."
          file_count = file_count + 1
        END DO
      END DO
    END DO
  END DO

  !$omp parallel private(kx, ky, k1, k2, band, Hamiltonian_const_band, Charge_dens_local, Delta_local,&
  !$omp&                 Gamma_K, orb, spin, n, neigh, lat, orb_prime, band_prime, layer, file_count, row, col, &
  !$omp&                 gamma_spin_index, gamma_lat_index, Hamiltonian_dummy, gamma_nn_12, gamma_nn_21, gamma_nnn)
  !$omp do
  DO i = -kx_steps, kx_steps
    DO j = -ky_steps, ky_steps
      kx = i * dkx
      ky = j * dky

      DO band = 1, SUBBANDS
        Hamiltonian_const_band = Hamiltonian_const
        CALL COMPUTE_SUBBAND_POTENTIAL(Hamiltonian_const_band, band)

        k1 = SQRT(3.0d0) / (2.0d0 * PI) * kx
        k2 = 3.0d0 / (4.0d0 * PI) * (ky + 1.0d0 / SQRT(3.0d0) * kx)
        CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const_band, Gamma_SC(:, :, :, :, band), Charge_dens(:, band), k1, k2, Delta_local(:, :, :, :, band), Charge_dens_local(:, band))
      END DO

      CALL GET_GAMMAS_FROM_DELTAS(Gamma_K, Delta_local)

      !Write result for given k point to file
      DO band = 1, SUBBANDS
        Hamiltonian_dummy = DCMPLX(0.0d0, 0.0d0)
        CALL COMPUTE_SC(Hamiltonian_dummy, kx, ky, Gamma_SC(:, :, :, :, band))

        DO orb = 1, ORBITALS
          DO spin = 0, 1
            gamma_spin_index = MOD(spin + 1, 2) + 1
            DO lat = 0, SUBLATTICES - 2
              !Include inter and intralayer couplings

              !Interlayer coupling
              !Coupling Ti1 - Ti2
              gamma_lat_index = 2 * lat + 2
              row = spin * TBA_DIM + orb + lat * ORBITALS
              col = orb + (lat + 1) * ORBITALS + DIM_POSITIVE_K + TBA_DIM * MOD(spin + 1, 2)
              gamma_nn_12 = Hamiltonian_dummy(row, col)

              IF (gamma_lat_index .le. SUBLATTICES) THEN
                row = orb + (gamma_lat_index - 1) * ORBITALS + spin * TBA_DIM
                col = orb + (gamma_lat_index - 1) * ORBITALS + MOD(spin + 1, 2) * TBA_DIM + DIM_POSITIVE_K
                gamma_nnn = Hamiltonian_dummy(row, col)
              ELSE
                gamma_nnn = DCMPLX(0.0d0, 0.0d0)
              END IF

              file_count = File_unit_mapping(orb, gamma_spin_index, gamma_lat_index, band)
              WRITE (file_count, '(6F15.8, *(2F15.8))') kx, ky, &
              & REAL(gamma_nn_12) / meV2au, AIMAG(gamma_nn_12) / meV2au, & !Interlayer couplings
              & REAL(gamma_nnn) / meV2au, AIMAG(gamma_nnn) / meV2au, & !Intralayer couplings
              & (REAL(Gamma_K(orb, neigh, gamma_spin_index, gamma_lat_index, band)) / meV2au, AIMAG(Gamma_K(orb, neigh, gamma_spin_index, gamma_lat_index, band)) / meV2au, neigh=1, N_ALL_NEIGHBOURS)

              !Coupling Ti2 - Ti1
              gamma_lat_index = 2 * lat + 1
              row = spin * TBA_DIM + orb + (lat + 1) * ORBITALS
              col = orb + DIM_POSITIVE_K + TBA_DIM * MOD(spin + 1, 2) + lat * ORBITALS
              gamma_nn_21 = Hamiltonian_dummy(row, col)

              IF (gamma_lat_index .le. SUBLATTICES) THEN
                row = orb + (gamma_lat_index - 1) * ORBITALS + spin * TBA_DIM
                col = orb + (gamma_lat_index - 1) * ORBITALS + MOD(spin + 1, 2) * TBA_DIM + DIM_POSITIVE_K
                gamma_nnn = Hamiltonian_dummy(row, col)
              ELSE
                gamma_nnn = DCMPLX(0.0d0, 0.0d0)
              END IF

              file_count = File_unit_mapping(orb, gamma_spin_index, gamma_lat_index, band)
              WRITE (file_count, '(6F15.8, *(2F15.8))') kx, ky, &
              & REAL(gamma_nn_21) / meV2au, AIMAG(gamma_nn_21) / meV2au, & !Interlayer couplings
              & REAL(gamma_nnn) / meV2au, AIMAG(gamma_nnn) / meV2au, & !Intralayer couplings
              & (REAL(Gamma_K(orb, neigh, gamma_spin_index, gamma_lat_index, band)) / meV2au, AIMAG(Gamma_K(orb, neigh, gamma_spin_index, gamma_lat_index, band)) / meV2au, neigh=1, N_ALL_NEIGHBOURS)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !$omp end do
  !$omp end parallel

  !Close all files
  DO orb = 1, ORBITALS
    DO spin = 1, 2
      DO layer = 1, LAYER_COUPLINGS
        DO band = 1, SUBBANDS
          file_count = File_unit_mapping(orb, spin, layer, band)
          CLOSE (file_count)
        END DO
      END DO
    END DO
  END DO

  DEALLOCATE (Hamiltonian_const)
  DEALLOCATE (Hamiltonian_const_band)
  DEALLOCATE (Hamiltonian_dummy)
  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Gamma_K)
  DEALLOCATE (Charge_dens)
  DEALLOCATE (Delta_local)
  DEALLOCATE (Charge_dens_local)
  DEALLOCATE (File_unit_mapping)
END SUBROUTINE CALCULATE_GAMMA_K

SUBROUTINE CALCULATE_PROJECTIONS(input_path, n_r_points, n_phi_points)
  CHARACTER(LEN=*), INTENT(IN) :: input_path !! This should be a path to folder where input.nml resides
  INTEGER, INTENT(IN) :: n_r_points !! Number of steps in radial direction
  INTEGER, INTENT(IN) :: n_phi_points !! Number of steps taken in angular direction in one of six triangles that make up the hexagon

  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :, :) !! Supercondicting pairings read from a simulation
  REAL*8, ALLOCATABLE :: Charge_dens(:, :) !! Charge density read from a simulation

  REAL*8 :: Kappa_nearest(3) !! K-variables aligned with orbital's directions for nearest neighbours
  REAL*8 :: Kappa_next(3) !! K-variables aligned with orbital's directions for next nearest neighbours
  COMPLEX*16 :: C_l(3) !! Shape functions for each orbital
  COMPLEX*16 :: Gamma_nearest_orb(3) !! Fourier transformed Gamma for given orbital at current k-point for nearest neighbour
  COMPLEX*16 :: Gamma_next_orb(3) !! Fourier transformed Gamma for given orbital at current k-point for next nearest neighbour
  INTEGER*8, PARAMETER :: N_PROJECTIONS = 18 !! Number of distict projections in orbital \otimes spatial basis
  COMPLEX*16 :: Projections_k_space_nearest(N_PROJECTIONS) = DCMPLX(0.0d0, 0.0d0) !! Accumulator for integrated projections over k-space for nearest neighbours
  COMPLEX*16 :: Projections_k_space_next(N_PROJECTIONS) = DCMPLX(0.0d0, 0.0d0) !! Accumulator for integrated projections over k-space for next nearest neighbours
  COMPLEX*16 :: Projections_real_space(N_PROJECTIONS) = DCMPLX(0.0d0, 0.0d0) !! Projections onto irreducible representations from real-space Gammas
  COMPLEX*16 :: Gamma_flat(N_PROJECTIONS)!! Flattened array of gammas in real space, to be used in real-space resolved symmetries
  COMPLEX*16 :: Active_orbital(3) !! Orbital for which to write projection. Used to utilize projection callbacks in writing basis functions for each orbital

  !File operations
  CHARACTER(LEN=200) :: filename
  INTEGER*4, PARAMETER:: FIRST_FILE_UNIT = 10
  INTEGER*4 :: file_count, gamma_nearest_weighted_file, gamma_next_weighted_file
  CHARACTER(LEN=4) :: Projections_name_mapping(N_PROJECTIONS)

  !K-space variables
  INTEGER*4 :: n_triangle, i_r, j_phi
  REAL*8 :: phi_k, r_k, r_max, dr
  REAL*8 :: kx, ky

  !Iterators
  INTEGER*4 :: n, neigh, orb, spin, layer, band

  ! Interface for basis functions
  ABSTRACT INTERFACE
    FUNCTION projection_interface(K_orb, Gamma_projected) RESULT(projection)
      REAL*8, INTENT(IN) :: K_orb(3)
      COMPLEX*16, INTENT(IN) :: Gamma_projected(3)
      COMPLEX*16 :: projection
    END FUNCTION projection_interface
  END INTERFACE
  ! Define a pointer to function type
  TYPE projection_function_cb_t
    PROCEDURE(projection_interface), POINTER, NOPASS :: cb
  END TYPE projection_function_cb_t
  ! Array of basis function pointers
  TYPE(projection_function_cb_t) :: Projections_cb(N_PROJECTIONS)

  CALL ASSIGN_PROJECTION_CALLBACKS(Projections_cb)
  CALL ASSIGN_PROJECTION_NAMES(Projections_name_mapping)

  CALL GET_INPUT(TRIM(input_path)//"input.nml")

  !Redefining step in radial and angular directions according to postprocessing parameters
  dr_k = R_K_MAX / n_r_points
  dphi_k = (PI / 3.0d0) / n_phi_points  !Slicing every hexagon's triangle into the same number of phi steps

  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K, SUBBANDS))

  CALL GET_SAFE_GAMMA_SC(Gamma_SC, input_path)
  CALL GET_SAFE_CHARGE_DENS(Charge_dens, input_path)

  Gamma_SC = 0.0d0
  Gamma_SC(:,:,1,:,:) = 1*meV2au
  Gamma_SC(:,:,2,:,:) = -1*meV2au

  CALL GET_REAL_SPACE_PROJECTIONS(Projections_real_space, Gamma_SC)
  DO n = 1, N_PROJECTIONS
    PRINT *, "Projection real space ", TRIM(ADJUSTL(Projections_name_mapping(n))), " = ", Projections_real_space(n) / SUM(Projections_real_space)
  END DO

  !Open files for projections
  file_count = FIRST_FILE_UNIT
  !Open files for projection basis functions
  DO n = 1, N_PROJECTIONS
    WRITE (filename, "(A, A)") "BasisFunc_", TRIM(ADJUSTL(Projections_name_mapping(n)))
    OPEN (unit=file_count, FILE=TRIM(input_path)//"OutputData/"//TRIM(filename)//".dat", FORM="FORMATTED", ACTION="WRITE")
    WRITE (file_count, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_basis_orb_nearest_1)[meV]   Im(Gamma_basis_orb_nearest_1)[meV]   Re(Gamma_basis_orb_next_1)[meV]   Im(Gamma_basis_orb_next_1)[meV]  Re(Gamma_basis_orb_nearest_2)[meV] ... Re(Gamma_basis_weighted)[meV]   Im(Gamma_basis_weighted)[meV]"
    file_count = file_count + 1
  END DO
  gamma_nearest_weighted_file = file_count
  file_count = file_count + 1
  OPEN (unit=gamma_nearest_weighted_file, FILE=TRIM(input_path)//"OutputData/Gamma_K_nearest_weighted.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_nearest_weighted_file, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_orb_1)[meV]   Im(Gamma_orb_1)[meV]  Re(Gamma_basis_2)[meV] ... Re(Gamma_weighted)[meV]   Im(Gamma_weighted)[meV]"

  gamma_next_weighted_file = file_count
  OPEN (unit=gamma_next_weighted_file, FILE=TRIM(input_path)//"OutputData/Gamma_K_next_weighted.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_next_weighted_file, '(A)') "#kx[1/a]   ky[1/a]   Re(Gamma_orb_1)[meV]   Im(Gamma_orb_1)[meV]  Re(Gamma_basis_2)[meV] ... Re(Gamma_weighted)[meV]   Im(Gamma_weighted)[meV]"

  !$omp parallel do collapse(3) schedule(dynamic, 1) private(phi_k, r_k, r_max, dr, kx, ky, Kappa_nearest, Kappa_next, C_l, Gamma_nearest_orb, Gamma_next_orb, n_triangle, j_phi, i_r, orb, n, Active_orbital)
  DO n_triangle = -N_BZ_SECTIONS / 2, N_BZ_SECTIONS / 2 - 1
    DO j_phi = 0, n_phi_points - 1
      DO i_r = 0, n_r_points - 1
        !Must stuck everhing here to provide pure loops for collapse(3)
        phi_k = n_triangle * (PI / 3.0d0) + j_phi * dphi_k
        r_max = r_max_phi(MOD(ABS(phi_k), PI / 3))
        dr = r_max / n_r_points
        r_k = i_r * dr

        !Transform from graphene reciprocal lattice to kx and ky
        kx = r_k * COS(phi_k)
        ky = r_k * SIN(phi_k)

        !Calculate k-space coordinates aligned with orbitals' directions
        Kappa_nearest(1) = SQRT(3.0d0) / 2.0d0 * kx + ky / 2.0d0
        Kappa_nearest(2) = -SQRT(3.0d0) / 2.0d0 * kx + ky / 2.0d0
        Kappa_nearest(3) = -ky

        Kappa_next(1) = -SQRT(3.0d0) / 2.0d0 * kx + 3.0d0 * ky / 2.0d0
        Kappa_next(2) = -SQRT(3.0d0) / 2.0d0 * kx - 3.0d0 * ky / 2.0d0
        Kappa_next(3) = SQRT(3.0d0) * kx

        !Calculate shape functions for orbitals
        C_l(1) = (1 + Kappa_nearest(2) * Kappa_nearest(3)) / (1 + kx**2 + ky**2)
        C_l(2) = (1 + Kappa_nearest(1) * Kappa_nearest(3)) / (1 + kx**2 + ky**2)
        C_l(3) = (1 + Kappa_nearest(1) * Kappa_nearest(2)) / (1 + kx**2 + ky**2)

        Gamma_nearest_orb = DCMPLX(0.0d0, 0.0d0)
        Gamma_next_orb = DCMPLX(0.0d0, 0.0d0)
        !Performing simple Fourier transform of Gamma in each orbital
        DO orb = 1, ORBITALS
          !Nearest neighbours
          Gamma_nearest_orb(orb) = Gamma_SC(orb, 1, 1, 1, 1) * pairing_1(ky) + Gamma_SC(orb, 1, 1, 2, 1) * CONJG(pairing_1(ky)) + &
                                   Gamma_SC(orb, 2, 1, 1, 1) * pairing_2(kx, ky) + Gamma_SC(orb, 2, 1, 2, 1) * CONJG(pairing_2(kx, ky)) + &
                                   Gamma_SC(orb, 3, 1, 1, 1) * pairing_3(kx, ky) + Gamma_SC(orb, 3, 1, 2, 1) * CONJG(pairing_3(kx, ky))

          !Next-nearest neighbours
          Gamma_next_orb(orb) = Gamma_SC(orb, N_NEIGHBOURS + 1, 1, 1, 1) * pairing_nnn_1(kx) + Gamma_SC(orb, N_NEIGHBOURS + 2, 1, 1, 1) * pairing_nnn_2(kx, ky) + &
                                Gamma_SC(orb, N_NEIGHBOURS + 3, 1, 1, 1) * pairing_nnn_3(kx, ky) + Gamma_SC(orb, N_NEIGHBOURS + 4, 1, 1, 1) * pairing_nnn_4(kx) + &
                                Gamma_SC(orb, N_NEIGHBOURS + 5, 1, 1, 1) * pairing_nnn_5(kx, ky) + Gamma_SC(orb, N_NEIGHBOURS + 6, 1, 1, 1) * pairing_nnn_6(kx, ky)
        END DO

        !Write Gamma(k) resulting from the file that was read
        WRITE (gamma_nearest_weighted_file, '(10F15.5)') kx, ky, (REAL(Gamma_nearest_orb(orb)), AIMAG(Gamma_nearest_orb(orb)), orb=1, ORBITALS),&
        & REAL(SUM(C_l * Gamma_nearest_orb)), AIMAG(SUM(C_l * Gamma_nearest_orb))

        WRITE (gamma_next_weighted_file, '(10F15.5)') kx, ky, (REAL(Gamma_next_orb(orb)), AIMAG(Gamma_next_orb(orb)), orb=1, ORBITALS),&
        & REAL(SUM(C_l * Gamma_next_orb)), AIMAG(SUM(C_l * Gamma_next_orb))

        !Project onto basis functions and write basis functions onto which we project
        DO n = 1, N_PROJECTIONS
          !$omp critical (projections_update)
          Projections_k_space_nearest(n) = Projections_k_space_nearest(n) + Projections_cb(n) % cb(Kappa_nearest, Gamma_nearest_orb) * r_k * dr * dphi_k
          Projections_k_space_next(n) = Projections_k_space_next(n) + Projections_cb(n) % cb(Kappa_next, Gamma_next_orb) * r_k * dr * dphi_k
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') kx, ky
          DO orb = 1, ORBITALS
            Active_orbital = 0.0d0
            Active_orbital(orb) = 1.0d0
            WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_nearest, Active_orbital)), AIMAG(Projections_cb(n) % cb(Kappa_nearest, Active_orbital))
            WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_next, Active_orbital)), AIMAG(Projections_cb(n) % cb(Kappa_next, Active_orbital))
          END DO
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)', ADVANCE='NO') REAL(Projections_cb(n) % cb(Kappa_nearest, C_l)), AIMAG(Projections_cb(n) % cb(Kappa_nearest, C_l))
          WRITE (FIRST_FILE_UNIT + n - 1, '(2F15.5)') REAL(Projections_cb(n) % cb(Kappa_next, C_l)), AIMAG(Projections_cb(n) % cb(Kappa_next, C_l))
        END DO
        !$omp end critical (projections_update)
      END DO
    END DO
  END DO
  !$omp end parallel do

  DO n = 1, N_PROJECTIONS
    PRINT *, "Projection k-space ", TRIM(ADJUSTL(Projections_name_mapping(n))), " = ", Projections_k_space_nearest(n) / SUM(Projections_k_space_nearest)
  END DO

  !Writing normalized projections since we are interested only in relative contributions
  OPEN (unit=gamma_next_weighted_file + 1, FILE=TRIM(input_path)//"OutputData/Projections_nearest.dat", FORM="FORMATTED", ACTION="WRITE")
  WRITE (gamma_next_weighted_file + 1, '(A)') "#Irrep      ABS(R-space projection)      ABS(K-space projection)"
  DO n = 1, N_PROJECTIONS
    WRITE (gamma_next_weighted_file + 1, '(A, 2F15.5)') Projections_name_mapping(n), &
    & ABS(Projections_real_space(n) / SUM(Projections_real_space)), ABS(Projections_k_space_nearest(n) / SUM(Projections_k_space_nearest)), &
    & ABS(Projections_k_space_next(n) / SUM(Projections_k_space_next))
  END DO
  CLOSE (gamma_next_weighted_file + 1)

  CLOSE (gamma_nearest_weighted_file)
  CLOSE (gamma_next_weighted_file)
  DO n = 1, N_PROJECTIONS
    CLOSE (FIRST_FILE_UNIT + n - 1)
  END DO

  DEALLOCATE (Gamma_SC)
  DEALLOCATE (Charge_dens)

! Internal functions
CONTAINS

  SUBROUTINE ASSIGN_PROJECTION_CALLBACKS(Projections_callbacks)
    !! Assign callbacks to projection onto irreducible representations functions
    IMPLICIT NONE
    TYPE(projection_function_cb_t), INTENT(OUT) :: Projections_callbacks(N_PROJECTIONS) !! Array of projection callbacks

    !Assigning callbacks of projections
    Projections_callbacks(1) % cb => a1_1_proj
    Projections_callbacks(2) % cb => a1_2_proj
    Projections_callbacks(3) % cb => a2_1_proj
    Projections_callbacks(4) % cb => b1_1_proj
    Projections_callbacks(5) % cb => b1_2_proj
    Projections_callbacks(6) % cb => b2_1_proj
    Projections_callbacks(7) % cb => e1_1_proj
    Projections_callbacks(8) % cb => e1_2_proj
    Projections_callbacks(9) % cb => e1_3_proj
    Projections_callbacks(10) % cb => e1_4_proj
    Projections_callbacks(11) % cb => e1_5_proj
    Projections_callbacks(12) % cb => e1_6_proj
    Projections_callbacks(13) % cb => e2_1_proj
    Projections_callbacks(14) % cb => e2_2_proj
    Projections_callbacks(15) % cb => e2_3_proj
    Projections_callbacks(16) % cb => e2_4_proj
    Projections_callbacks(17) % cb => e2_5_proj
    Projections_callbacks(18) % cb => e2_6_proj
  END SUBROUTINE ASSIGN_PROJECTION_CALLBACKS

  SUBROUTINE ASSIGN_PROJECTION_NAMES(Projections_names)
    !! Map index of projection to its name
    IMPLICIT NONE
    CHARACTER(LEN=4) :: Projections_names(N_PROJECTIONS) !! Array of projection names mapped to projection number

    Projections_names(1) = "A1_1"
    Projections_names(2) = "A1_2"
    Projections_names(3) = "A2_1"
    Projections_names(4) = "B1_1"
    Projections_names(5) = "B1_2"
    Projections_names(6) = "B2_1"
    Projections_names(7) = "E1_1"
    Projections_names(8) = "E1_2"
    Projections_names(9) = "E1_3"
    Projections_names(10) = "E1_4"
    Projections_names(11) = "E1_5"
    Projections_names(12) = "E1_6"
    Projections_names(13) = "E2_1"
    Projections_names(14) = "E2_2"
    Projections_names(15) = "E2_3"
    Projections_names(16) = "E2_4"
    Projections_names(17) = "E2_5"
    Projections_names(18) = "E2_6"

  END SUBROUTINE ASSIGN_PROJECTION_NAMES

  SUBROUTINE GET_REAL_SPACE_PROJECTIONS(Projections_real_space, Gamma_SC)
    IMPLICIT NONE
    COMPLEX*16, INTENT(OUT) :: Projections_real_space(N_PROJECTIONS) !! Array of projections in real space
    COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS, SUBBANDS)
    COMPLEX*16 :: Gamma_flat(N_PROJECTIONS)
    INTEGER*4 :: orb, neigh, n
    n = 1
    DO orb = 1, ORBITALS
      DO neigh = 1, N_NEIGHBOURS
        Gamma_flat(n) = Gamma_SC(orb, neigh, 1, 1, 1)
        n = n + 1
        !Get neighbor from opposite sublattice corresponding to the next vector in conter-clockwise rotation
        Gamma_flat(n) = Gamma_SC(orb, MOD(neigh + 1, N_NEIGHBOURS) + 1, 1, 2, 1)
        n = n + 1
      END DO
    END DO
    Projections_real_space(1) = (Gamma_flat(3) + Gamma_flat(6) + Gamma_flat(8) + Gamma_flat(11) + Gamma_flat(13) + Gamma_flat(16))
    Projections_real_space(2) = (Gamma_flat(1) + Gamma_flat(2) + Gamma_flat(4) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(9) + Gamma_flat(10) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(15) + Gamma_flat(17) + Gamma_flat(18))
    Projections_real_space(3) = (Gamma_flat(2) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(10) + Gamma_flat(15) + Gamma_flat(18) - (Gamma_flat(1) + Gamma_flat(4) + Gamma_flat(9) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(17)))
    Projections_real_space(4) = (Gamma_flat(6) + Gamma_flat(8) + Gamma_flat(16) - (Gamma_flat(3) + Gamma_flat(11) + Gamma_flat(13)))
    Projections_real_space(5) = (Gamma_flat(2) + Gamma_flat(4) + Gamma_flat(10) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(18) - (Gamma_flat(1) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(9) + Gamma_flat(15) + Gamma_flat(17)))
    Projections_real_space(6) = (Gamma_flat(1) + Gamma_flat(2) + Gamma_flat(9) + Gamma_flat(10) + Gamma_flat(17) + Gamma_flat(18) - (Gamma_flat(4) + Gamma_flat(5) + Gamma_flat(7) + Gamma_flat(12) + Gamma_flat(14) + Gamma_flat(15)))
    Projections_real_space(7) = (Gamma_flat(5) + Gamma_flat(10) - (Gamma_flat(2) + Gamma_flat(7)))
    Projections_real_space(8) = (Gamma_flat(6) + Gamma_flat(11) - (Gamma_flat(3) + Gamma_flat(8)))
    Projections_real_space(9) = (Gamma_flat(1) + Gamma_flat(12) - (Gamma_flat(4) + Gamma_flat(9)))
    Projections_real_space(10) = (Gamma_flat(3) + Gamma_flat(16) - (Gamma_flat(6) + Gamma_flat(13)))
    Projections_real_space(11) = (Gamma_flat(4) + Gamma_flat(17) - (Gamma_flat(1) + Gamma_flat(14)))
    Projections_real_space(12) = (Gamma_flat(5) + Gamma_flat(18) - (Gamma_flat(2) + Gamma_flat(15)))
    Projections_real_space(13) = (Gamma_flat(7) + Gamma_flat(10) - Gamma_flat(2) - Gamma_flat(5))
    Projections_real_space(14) = (Gamma_flat(8) + Gamma_flat(11) - Gamma_flat(3) - Gamma_flat(6))
    Projections_real_space(15) = (Gamma_flat(9) + Gamma_flat(12) - Gamma_flat(1) - Gamma_flat(4))
    Projections_real_space(16) = (Gamma_flat(13) + Gamma_flat(16) - Gamma_flat(3) - Gamma_flat(6))
    Projections_real_space(17) = (Gamma_flat(14) + Gamma_flat(17) - Gamma_flat(1) - Gamma_flat(4))
    Projections_real_space(18) = (Gamma_flat(15) + Gamma_flat(18) - Gamma_flat(2) - Gamma_flat(5))

    !Normalize to per-bond-per-orbital coupling in meV
    Projections_real_space = Projections_real_space / N_PROJECTIONS

  END SUBROUTINE GET_REAL_SPACE_PROJECTIONS

  !dir$ attributes forceinline :: a1_1_proj
  PURE COMPLEX * 16 FUNCTION a1_1_proj(k_orb, Gamma_projected)
    !! Project Gamma_projected onto A1_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(3) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * COS(K_orb(1))
    F_orb(2) = 2 * COS(K_orb(2))
    F_orb(3) = 2 * COS(K_orb(3))
    a1_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION a1_1_proj

  !dir$ attributes forceinline :: a1_2_proj
  PURE COMPLEX * 16 FUNCTION a1_2_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto A1_2 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS)
    F_orb(1) = 2 * (COS(K_orb(2)) + COS(K_orb(3)))
    F_orb(2) = 2 * (COS(K_orb(1)) + COS(K_orb(3)))
    F_orb(3) = 4 * COS(K_orb(3) / 2.0d0) * COS((K_orb(1) - K_orb(2)) / 2.0d0)
    a1_2_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION a1_2_proj

  !dir$ attributes forceinline :: a2_1_proj
  PURE COMPLEX * 16 FUNCTION a2_1_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto A2_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * (COS(K_orb(2)) - COS(K_orb(3)))
    F_orb(2) = 2 * (-COS(K_orb(1)) + COS(K_orb(3)))
    F_orb(3) = 4 * SIN(K_orb(3) / 2.0d0) * SIN((K_orb(1) - K_orb(2)) / 2.0d0)
    a2_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION a2_1_proj

  !dir$ attributes forceinline :: b1_1_proj
  PURE COMPLEX * 16 FUNCTION b1_1_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto B1_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(1))
    F_orb(2) = 2 * imag * SIN(K_orb(2))
    F_orb(3) = 2 * imag * SIN(K_orb(3))
    b1_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION b1_1_proj

  !dir$ attributes forceinline :: b1_2_proj
  PURE COMPLEX * 16 FUNCTION b1_2_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto B1_2 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * (SIN(K_orb(2)) + SIN(K_orb(3)))
    F_orb(2) = 2 * imag * (SIN(K_orb(1)) + SIN(K_orb(3)))
    F_orb(3) = -4 * imag * SIN(K_orb(3) / 2.0d0) * COS((K_orb(1) - K_orb(2)) / 2.0d0)
    b1_2_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION b1_2_proj

  !dir$ attributes forceinline :: b2_1_proj
  PURE COMPLEX * 16 FUNCTION b2_1_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto B2_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * (SIN(K_orb(2)) - SIN(K_orb(3)))
    F_orb(2) = 2 * imag * (-SIN(K_orb(1)) + SIN(K_orb(3)))
    F_orb(3) = 4 * imag * COS(K_orb(3) / 2.0d0) * SIN((K_orb(1) - K_orb(2)) / 2.0d0)
    b2_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION b2_1_proj

  !dir$ attributes forceinline :: e1_1_proj
  PURE COMPLEX * 16 FUNCTION e1_1_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(2))
    F_orb(2) = 2 * imag * SIN(K_orb(3))
    F_orb(3) = 0.d0
    e1_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_1_proj

  !dir$ attributes forceinline :: e1_2_proj
  PURE COMPLEX * 16 FUNCTION e1_2_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_2 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(1))
    F_orb(2) = -2 * imag * SIN(K_orb(2))
    F_orb(3) = 0.d0
    e1_2_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_2_proj

  !dir$ attributes forceinline :: e1_3_proj
  PURE COMPLEX * 16 FUNCTION e1_3_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_3 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(3))
    F_orb(2) = 2 * imag * SIN(K_orb(1))
    F_orb(3) = 0.d0
    e1_3_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_3_proj

  !dir$ attributes forceinline :: e1_4_proj
  PURE COMPLEX * 16 FUNCTION e1_4_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_4 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(1))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * imag * SIN(K_orb(3))
    e1_4_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_4_proj

  !dir$ attributes forceinline :: e1_5_proj
  PURE COMPLEX * 16 FUNCTION e1_5_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_5 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = 2 * imag * SIN(K_orb(3))
    F_orb(2) = 0.d0
    F_orb(3) = -2 * imag * SIN(K_orb(2))
    e1_5_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_5_proj

  !dir$ attributes forceinline :: e1_6_proj
  PURE COMPLEX * 16 FUNCTION e1_6_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E1_6 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * imag * SIN(K_orb(2))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * imag * SIN(K_orb(1))
    e1_6_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e1_6_proj

  !dir$ attributes forceinline :: e2_1_proj
  PURE COMPLEX * 16 FUNCTION e2_1_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_1 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(2))
    F_orb(2) = 2 * COS(K_orb(3))
    F_orb(3) = 0.d0
    e2_1_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_1_proj

  !dir$ attributes forceinline :: e2_2_proj
  PURE COMPLEX * 16 FUNCTION e2_2_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_2 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(1))
    F_orb(2) = 2 * COS(K_orb(2))
    F_orb(3) = 0.d0
    e2_2_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_2_proj

  !dir$ attributes forceinline :: e2_3_proj
  PURE COMPLEX * 16 FUNCTION e2_3_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_3 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(3))
    F_orb(2) = 2 * COS(K_orb(1))
    F_orb(3) = 0.d0
    e2_3_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_3_proj

  !dir$ attributes forceinline :: e2_4_proj
  PURE COMPLEX * 16 FUNCTION e2_4_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_4 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(1))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(3))
    e2_4_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_4_proj

  !dir$ attributes forceinline :: e2_5_proj
  PURE COMPLEX * 16 FUNCTION e2_5_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_5 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(3))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(2))
    e2_5_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_5_proj

  !dir$ attributes forceinline :: e2_6_proj
  PURE COMPLEX * 16 FUNCTION e2_6_proj(K_orb, Gamma_projected)
    !! Project Gamma_projected onto E2_6 irreducible representation
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: K_orb(ORBITALS) !! Set of orbital-aligned k-space coordinates
    COMPLEX*16, INTENT(IN) :: Gamma_projected(ORBITALS) !! Set of k-dependent Gammas to be projected
    COMPLEX*16 :: F_orb(ORBITALS) !! Basis function for each orbital at given irreducible representation
    F_orb(1) = -2 * COS(K_orb(2))
    F_orb(2) = 0.d0
    F_orb(3) = 2 * COS(K_orb(1))
    e2_6_proj = DOT_PRODUCT(F_orb, Gamma_projected)
    RETURN
  END FUNCTION e2_6_proj

END SUBROUTINE CALCULATE_PROJECTIONS

SUBROUTINE TRANSFORM_DELTA_MATRIX(inputPath, nBrillouinPoints)
  CHARACTER(LEN=*), INTENT(IN) :: inputPath !! This should be a path to folder where input.nml resides
  INTEGER*4, INTENT(IN) :: nBrillouinPoints

  CHARACTER(LEN=20) :: output_format

  COMPLEX*16, ALLOCATABLE :: Hamiltonian(:, :), Hamiltonian_const(:, :), Gamma_matrix(:, :), Gamma_matrix_temp(:, :), Gamma_matrix_diagonal(:, :)
  COMPLEX*16, ALLOCATABLE :: U_transformation(:, :)
  REAL*8, ALLOCATABLE :: Energies(:)

  COMPLEX*16, ALLOCATABLE :: Gamma_SC(:, :, :, :)
  REAL*8, ALLOCATABLE :: Charge_dens(:)
  REAL*8 :: brillouinZoneVertices(6, 2)

  REAL*8 :: k1, k2, kx, ky, dkx, dky
  INTEGER*4 :: i, j, n, m
  INTEGER*4 :: kx_steps, ky_steps

  LOGICAL :: fileExists

  brillouinZoneVertices(:, 1) = (/4.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), -4.*PI / (3 * SQRT(3.0d0)), -2.*PI / (3 * SQRT(3.0d0)), 2.*PI / (3 * SQRT(3.0d0))/)
  brillouinZoneVertices(:, 2) = (/0.0d0, -2.*PI / 3.0d0, -2.*PI / 3.0d0, 0.0d0, 2.*PI / 3.0d0, 2.*PI / 3.0d0/)

  dkx = KX_MAX / nBrillouinPoints
  dky = KY_MAX / nBrillouinPoints

  kx_steps = INT(nBrillouinPoints)
  ky_steps = INT(nBrillouinPoints)

  ALLOCATE (Hamiltonian(DIM, DIM))
  ALLOCATE (Hamiltonian_const(DIM, DIM))
  ALLOCATE (U_transformation(DIM, DIM))
  ALLOCATE (Gamma_matrix(DIM, DIM))
  ALLOCATE (Gamma_matrix_diagonal(DIM, DIM))
  ALLOCATE (Gamma_matrix_temp(DIM, DIM))
  ALLOCATE (Energies(DIM))
  ALLOCATE (Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, SUBLATTICES))
  ALLOCATE (Charge_dens(DIM_POSITIVE_K))

  Gamma_SC = DCMPLX(0.0d0, 0.0d0)
  Charge_dens = 0.0d0

  !Get parameters from simulation
  CALL GET_INPUT(TRIM(inputPath)//"input.nml")
  !Charge densities
  ! INQUIRE(FILE = TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST = fileExists)
  ! IF (fileExists) THEN
  !     CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
  ! ELSE
  !     CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
  ! END IF
  ! !Superconducting gap parameters
  ! INQUIRE(FILE = TRIM(inputPath)//"OutputData/Gamma_SC__final.dat", EXIST = fileExists)
  ! IF (fileExists) THEN
  !     CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
  ! ELSE
  !     CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
  ! END IF

  Gamma_SC(:, :, 1, :) = 999.*meV2au
  Gamma_SC(:, :, 2, :) = -999.*meV2au

  !Initialize constant hamiltonian
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  DO i = -kx_steps, kx_steps
    DO j = -ky_steps, ky_steps
      kx = i * dkx !* (2. * PI * 2./3.)
      ky = j * dky !* (2. * PI * 2./3.)

      IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

        Hamiltonian(:, :) = DCMPLX(0., 0.)
        Energies(:) = 0.
        CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky)
        CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:))
        CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :))
        CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

        Hamiltonian(:, :) = 0.5 * (Hamiltonian_const(:, :) + Hamiltonian(:, :)) !Should by multiplied by 0.5 if in Nambu space
        ! DO n = 1, DIM
        !     WRITE(*,'(24E12.4)') (ABS(Hamiltonian(n,m))**2/meV2au, m = 1, DIM)
        ! END DO
        ! WRITE(*,*)
        ! WRITE(*,*)

        !CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:, :), Energies(:), DIM)
        CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:, :), Energies(:), U_transformation(:, :), DIM)

        !Compute current superconducing coupling matrix separately
        Gamma_matrix(:, :) = DCMPLX(0.0d0, 0.0d0)
        Gamma_matrix_diagonal(:, :) = DCMPLX(0.0d0, 0.0d0)
        Gamma_matrix_temp(:, :) = DCMPLX(0.0d0, 0.0d0)
        CALL COMPUTE_SC(Gamma_matrix(:, :), kx, ky, Gamma_SC(:, :, :, :))
        CALL COMPUTE_CONJUGATE_ELEMENTS(Gamma_matrix(:, :), DIM)
        Gamma_matrix = 0.5 * Gamma_matrix
        DO n = 1, DIM
          DO m = 1, DIM
            Gamma_matrix_temp(n, m) = SUM(Gamma_matrix(n, :) * U_transformation(:, m))
          END DO
        END DO

        DO n = 1, DIM
          DO m = 1, DIM
            Gamma_matrix_diagonal(n, m) = SUM(CONJG(U_transformation(:, n)) * Gamma_matrix_temp(:, m))
          END DO
        END DO

        !Gamma_matrix_diagonal = MATMUL(MATMUL(Gamma_matrix, U_transformation), TRANSPOSE(CONJG(U_transformation)))

        WRITE (*, '(A, 24F10.4)') 'Energies: ', (Energies(m) / meV2au, m=1, DIM)
        DO n = 1, DIM
          WRITE (*, '(24F10.4)') (REAL(Gamma_matrix_diagonal(n, m)) / meV2au, m=1, DIM)
        END DO
        WRITE (*, *)
        WRITE (*, *)

      END IF
    END DO
  END DO

END SUBROUTINE TRANSFORM_DELTA_MATRIX

SUBROUTINE HELLICAL_TEST_CHERN(potChem, B, Nk1, Nk2, i, j, U_transformation)
  INTEGER*4, PARAMETER :: HamDim = 4
  REAL*8, PARAMETER :: g = 5.0d0
  REAL*8, PARAMETER :: tHop = 200 * meV2au
  REAL*8, PARAMETER :: alphaSOC = 100 * meV2au
  REAL*8, PARAMETER :: gammaSC = DCMPLX(0.5 * meV2au, 0.0d0)
  REAL*8, PARAMETER :: muB = 0.5d0

  REAL*8, INTENT(IN) :: potChem
  REAL*8, INTENT(IN) :: B(3) ![Bx, By, Bz]
  INTEGER*4, INTENT(IN) :: Nk1, Nk2, i, j

  COMPLEX*16, INTENT(INOUT) :: U_transformation(HamDim, HamDim)
  REAL*8 :: Energies(HamDim)
  REAL*8 :: dkx, dky
  REAL*8 :: kx, ky
  INTEGER*4 :: m, n

  COMPLEX*16 :: Hamiltonian(HamDim, HamDim)

  dkx = 2.0d0 * PI / Nk1
  dky = 2.0d0 * PI / Nk2

  Energies(:) = 0.0d0
  U_transformation(:, :) = DCMPLX(0.0d0, 0.0d0)
  Hamiltonian(:, :) = 0.0d0
  !Calculate all points in Brillouin zone
  kx = i * dkx
  ky = j * dky

  !Diagonal terms
  !Electrons H(k)
  Hamiltonian(1, 1) = 2.0 * tHop * (1.-DCOS(kx)) + 2.0 * tHop * (1 - DCOS(ky)) - potChem + 0.5 * muB * g * B(3)
  Hamiltonian(2, 2) = 2.0 * tHop * (1.-DCOS(kx)) + 2.0 * tHop * (1 - DCOS(ky)) - potChem - 0.5 * muB * g * B(3)

  !Holes -H*(-k)
  Hamiltonian(3, 3) = -(2.0 * tHop * (1.-DCOS(-kx)) + 2.0 * tHop * (1 - DCOS(-ky)) - potChem + 0.5 * muB * g * B(3))
  Hamiltonian(4, 4) = -(2.0 * tHop * (1.-DCOS(-kx)) + 2.0 * tHop * (1 - DCOS(-ky)) - potChem - 0.5 * muB * g * B(3))

  !Spin-orbit coupling
  Hamiltonian(1, 2) = 0.5 * mub * g * (B(1) - imag * B(2)) + alphaSOC * (DSIN(kx) + imag * DSIN(ky))
  Hamiltonian(3, 4) = -(0.5 * mub * g * (B(1) + imag * B(2)) + alphaSOC * (DSIN(-kx) - imag * DSIN(-ky)))

  !Superconductivity
  Hamiltonian(1, 4) = gammaSC
  Hamiltonian(2, 3) = -gammaSC

  Hamiltonian(:, :) = 0.5 * Hamiltonian(:, :)
  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, HamDim)

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:, :), Energies(:), U_transformation(:, :), HamDim)

  CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies(:), U_transformation(:, :), HamDim)

END SUBROUTINE HELLICAL_TEST_CHERN

SUBROUTINE LAO_STO_CHERN_ENERGIES(Nk1, Nk2, i, j, inputPath, U_transformation)
    !! This subroutine calculates energies and wavefunctions of LAO-STO in [111] direction.
    !! Returns sorted wavefunctions in (i,j) point of the Brillouin zone.
  INTEGER*4, INTENT(IN) :: Nk1 !! Number of divisions of Brillouin zone in direction k1.
  INTEGER*4, INTENT(IN) :: Nk2 !! Number of divisions of Brillouin zone in direction k2.
  INTEGER*4, INTENT(IN) :: i !! Curent point k1_i in the Brillouin zone
  INTEGER*4, INTENT(IN) :: j !! Curent point k2_j in the Brillouin zone
  CHARACTER(LEN=*), INTENT(IN) :: inputPath !! Directory of the run, where input.nml should be placed.
                                            !! It contains material information and physical parameter of calculation:
                                            !! Fermi energy, temperature etc.
  COMPLEX*16, INTENT(OUT) :: U_transformation(DIM, DIM) !! Matrix containing eigenvectors stored in consecutive columns.
                                                          !! On output sorted based on energies from lowest to highest.

  REAL*8 :: Energies(DIM) !! Eigenvalues of the hamiltonian

  COMPLEX*16 :: Hamiltonian(DIM, DIM), Hamiltonian_const(DIM, DIM)
  COMPLEX*16 :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, SUBLATTICES)
  REAL*8 :: Charge_dens(DIM_POSITIVE_K)
  REAL*8 :: k1, k2, kx, ky
  REAL*8 :: dk1_Chern, dk2_Chern
  INTEGER*4 :: n
  LOGICAL :: fileExists

  !PRINT*, "Allocation ended"
  dk1_Chern = K1_MAX / Nk1
  dk2_Chern = K2_MAX / Nk2

  U_transformation = DCMPLX(0., 0.)
  Energies = 0.
  Gamma_SC = 0.
  Charge_dens = 0.

  !PRINT*, "Entered chern energies"
  !Get parameters from simulation
  CALL GET_INPUT(TRIM(inputPath)//"input.nml")
  !CHeck if units are correct
  !Get Gammas from run
  INQUIRE (FILE=TRIM(inputPath)//"OutputData/Gamma_SC__final.dat", EXIST=fileExists)
  IF (fileExists) THEN
    CALL GET_GAMMA_SC(Gamma_SC(:, :, :, :), TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
  ELSE
    CALL GET_GAMMA_SC(Gamma_SC(:, :, :, :), TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
  END IF

  !Get cherge densities from file
  INQUIRE (FILE=TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST=fileExists)
  IF (fileExists) THEN
    CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
  ELSE
    CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
  END IF

  !Gamma_SC(:,:,1,:) = 100.0d0 * meV2au
  !Gamma_SC(:,:,2,:) = -100.0d0 * meV2au

  !Computing k-independent terms
  Hamiltonian_const = DCMPLX(0., 0.)
  CALL COMPUTE_K_INDEPENDENT_TERMS(Hamiltonian_const)

  !Calculate eigenvalues and eigenvectors to later compute chern numbers
  k1 = i * dk1_Chern
  k2 = j * dk2_Chern

  kx = 2.*PI / (SQRT(3.0d0)) * k1
  ky = -2.*PI / 3.*k1 + 4.*PI / 3.*k2
  Hamiltonian(:, :) = DCMPLX(0., 0.)
  CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky)
  CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:))
  CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :))

  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

  Hamiltonian(:, :) = 0.5 * (Hamiltonian_const(:, :) + Hamiltonian(:, :)) !Should by multiplied by 0.5 if in Nambu space

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:, :), Energies(:), U_transformation(:, :), DIM)

  CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, U_transformation, DIM)

  ! PRINT*, "ENERGIES FOR ", i, j, k1, k2
  ! DO n = 1, DIM
  !     PRINT*, Energies(n)
  ! END DO

END SUBROUTINE LAO_STO_CHERN_ENERGIES

!########################### HELPER FUNCTIONS ###########################
SUBROUTINE SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, Psi, HamDim)
  INTEGER*4, INTENT(IN) :: HamDim
  COMPLEX*16, INTENT(INOUT) :: Psi(HamDim, HamDim)
  REAL*8, INTENT(INOUT) :: Energies(HamDim)

  INTEGER*4 :: i, j
  REAL*8 :: tmpEnergy
  COMPLEX*16 :: tmpPsi(HamDim)

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

COMPLEX * 16 FUNCTION det(matrix, n)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: n
  COMPLEX*16, INTENT(IN) :: matrix(n, n)
  INTEGER*4 :: IPIV(n)
  INTEGER*4 :: info, i

  IPIV(:) = 0.0d0

  CALL ZGETRF(n, n, matrix, n, IPIV, info)
  CALL ZLAPMT(.TRUE., n, n, matrix, n, IPIV)
  det = 1.0d0
  DO i = 1, n
    det = det * matrix(i, i)
  END DO
  RETURN

END FUNCTION det

LOGICAL FUNCTION is_inside_polygon(verticesArray, nVertices, pointX, pointY)
    !! This function returns true if the point (pointX, pointY) is inside the polygon
    !! Where vertices are defined in array verticesArray, containing X and Y coordinates.
    !! The assumption is that the polygon is two-dimensional
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: nVertices !! Number of vertices of polygon
  REAL*8, INTENT(IN) :: verticesArray(nVertices, 2) !! X and Y coordinates of vertices
  REAL*8, INTENT(IN) :: pointX, pointY !! X and Y coordinates of point to be tested
  INTEGER*4 :: i, j
  REAL*8 :: xIntersection

  is_inside_polygon = .FALSE.

  DO i = 1, nVertices

    IF (pointX == verticesArray(i, 1) .AND. pointY == verticesArray(i, 2)) THEN
      is_inside_polygon = .TRUE.
      RETURN
    END IF

    j = MOD(i, nVertices) + 1 !! Index of next vertex
        !! Check if horizontal line of y = pointY can cross the line joining two vertices
    IF (pointY > MIN(verticesArray(i, 2), verticesArray(j, 2)) .AND. &
    & pointY <= MAX(verticesArray(i, 2), verticesArray(j, 2))) THEN
      IF (pointX <= MAX(verticesArray(i, 1), verticesArray(j, 1))) THEN
        IF (verticesArray(i, 2) /= verticesArray(j, 2)) THEN !To avoid division by zero
          !Calculates intersection point of line connecting two vertices and horizontal line y = pointY
          xIntersection = (pointY - verticesArray(i, 2)) * &
          & (verticesArray(j, 1) - verticesArray(i, 1)) / (verticesArray(j, 2) &
          & - verticesArray(i, 2)) + verticesArray(i, 1)

          IF (pointX <= xIntersection .OR. verticesArray(i, 1) == verticesArray(j, 1)) THEN
            is_inside_polygon = .NOT. is_inside_polygon
          END IF
        END IF
      END IF
    END IF
  END DO

  RETURN
END FUNCTION is_inside_polygon

!MAKE THIS A SUBROUTINE
!Plot this dispersion
! OPEN(unit = 9, FILE= "./OutputData/EnergiesHellical.dat", FORM = "FORMATTED", ACTION = "WRITE")
! DO n = -Nk1/2, Nk1/2
!     kx = n*dkx
!     ky = kx

!     ! Hamiltonian(1, 1) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) + 0.5*muB*g*B(3)
!     ! Hamiltonian(2, 2) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - 0.5*muB*g*B(3)

!     ! !Spin-orbit coupling
!     ! Hamiltonian(1,2) = 0.5d0*mub*g*(B(1) - imag*B(2)) + alphaSOC*(DSIN(kx) + imag*DSIN(ky))
!     ! CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:2, :2), 2)

!     ! CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:2, :2), Energies(:2), U_transformation(:2,:2), 2)

!     !Diagonal terms
!     !Electrons H(k)
!     Hamiltonian(1, 1) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem + 0.5*muB*g*B(3)
!     Hamiltonian(2, 2) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem - 0.5*muB*g*B(3)

!     !Holes -H*(-k)
!     Hamiltonian(3, 3) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem + 0.5*muB*g*B(3))
!     Hamiltonian(4, 4) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem - 0.5*muB*g*B(3))

!     !Spin-orbit coupling
!     Hamiltonian(1,2) = 0.5*mub*g*(B(1) - imag*B(2)) + alphaSOC*(DSIN(kx) + imag*DSIN(ky))
!     Hamiltonian(3,4) = -(0.5*mub*g*(B(1) + imag*B(2)) + alphaSOC*(DSIN(-kx) - imag*DSIN(-ky)))

!     !Superconductivity
!     Hamiltonian(1, 4) = gammaSC
!     Hamiltonian(2,3) = -gammaSC

!     Hamiltonian(:,:) = 0.5 * Hamiltonian(:,:)
!     CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, HamDim)

!     CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), HamDim)

!     DO m = 1, HamDim
!         WRITE(9,*) kx, Energies(m)/meV2au
!     END DO

! END DO
! CLOSE(9)

END MODULE mod_postprocessing
