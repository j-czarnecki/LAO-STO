MODULE mod_local_integrand
USE mod_parameters
USE mod_utilities
USE mod_hamiltonians
USE mod_writers
IMPLICIT NONE
CONTAINS

SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, Charge_dens, k1, k2, Delta_local, Charge_dens_local)
  COMPLEX*16, INTENT(IN) :: Hamiltonian_const(DIM, DIM)
  REAL*8, INTENT(IN) :: k1, k2
  COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K)

  COMPLEX*16, INTENT(OUT) :: Delta_local(ORBITALS, N_ALL_NEIGHBOURS, 2, LAYER_COUPLINGS)
  REAL*8, INTENT(OUT) :: Charge_dens_local(DIM_POSITIVE_K)

  COMPLEX*16 :: Hamiltonian(DIM, DIM)
  COMPLEX*16 :: U_transformation(DIM, DIM)
  REAL*8 :: Energies(DIM)
  REAL*8 :: kx, ky
  INTEGER*4 :: orb, lat, n, m, spin
  INTEGER*4 :: row, col
  REAL*8 :: occupation

  !Transform from graphene reciprocal lattice to kx and ky
  kx = k1 * COS(k2)
  ky = k1 * SIN(k2)

  Energies(:) = 0.
  Hamiltonian(:, :) = DCMPLX(0., 0.)
  U_transformation(:, :) = DCMPLX(0., 0.)
  CALL COMPUTE_K_DEPENDENT_TERMS(Hamiltonian, kx, ky)
  CALL COMPUTE_HUBBARD(Hamiltonian(:, :), Charge_dens(:))
  CALL COMPUTE_SC(Hamiltonian(:, :), kx, ky, Gamma_SC(:, :, :, :))

  CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:, :), DIM) !This is not needed, since ZHEEV takes only upper triangle

  Hamiltonian(:, :) = 0.5 * (Hamiltonian_const(:, :) + Hamiltonian(:, :))
  !U_transformation(:,:) = Hamiltonian(:,:)
  !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:), DIM)
  ! CALL PRINT_HAMILTONIAN(Hamiltonian(:,:))
  ! STOP 'Hamiltonian printed'

  CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:, :), Energies(:), U_transformation(:, :), DIM)
  !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U

  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Delta_local(:, :, :, :) = DCMPLX(0., 0.)
  !Self - consistent delta calculation
  DO orb = 1, ORBITALS
    !### NEAREST NEIGHBOURS PAIRING ###########################################
    !Electrons
    DO n = 1, DIM_POSITIVE_K
      occupation = fd_distribution(Energies(n), 0d0, T)
      DO spin = 0, 1
        DO lat = 1, LAYER_COUPLINGS, 2
          !Coupling from current to next layer
          row = orb + DIM_POSITIVE_K + spin * TBA_DIM + (lat / 2) * ORBITALS
          col = orb + TBA_DIM * MOD(spin + 1, 2) + ((lat + 1) / 2) * ORBITALS
          Delta_local(orb, 1, spin + 1, lat) = Delta_local(orb, 1, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_1(ky)
          Delta_local(orb, 2, spin + 1, lat) = Delta_local(orb, 2, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_2(kx, ky)
          Delta_local(orb, 3, spin + 1, lat) = Delta_local(orb, 3, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_3(kx, ky)

          !Coupling from next to current layer
          row = orb + DIM_POSITIVE_K + spin * TBA_DIM + ((lat + 1) / 2) * ORBITALS
          col = orb + TBA_DIM * MOD(spin + 1, 2) + (lat / 2) * ORBITALS
          Delta_local(orb, 1, spin + 1, lat + 1) = Delta_local(orb, 1, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_1(ky))
          Delta_local(orb, 2, spin + 1, lat + 1) = Delta_local(orb, 2, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_2(kx, ky))
          Delta_local(orb, 3, spin + 1, lat + 1) = Delta_local(orb, 3, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_3(kx, ky))
        END DO
      END DO
    END DO

    !Holes
    DO n = DIM_POSITIVE_K + 1, DIM
      occupation = 1.-fd_distribution(-Energies(n), 0d0, T)
      DO spin = 0, 1
        DO lat = 1, LAYER_COUPLINGS, 2
          !Coupling from current to next layer
          row = orb + DIM_POSITIVE_K + spin * TBA_DIM + (lat / 2) * ORBITALS
          col = orb + TBA_DIM * MOD(spin + 1, 2) + ((lat + 1) / 2) * ORBITALS
          Delta_local(orb, 1, spin + 1, lat) = Delta_local(orb, 1, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_1(ky)
          Delta_local(orb, 2, spin + 1, lat) = Delta_local(orb, 2, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_2(kx, ky)
          Delta_local(orb, 3, spin + 1, lat) = Delta_local(orb, 3, spin + 1, lat) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * pairing_3(kx, ky)

          !Coupling from next to current layer
          row = orb + DIM_POSITIVE_K + spin * TBA_DIM + ((lat + 1) / 2) * ORBITALS
          col = orb + TBA_DIM * MOD(spin + 1, 2) + (lat / 2) * ORBITALS
          Delta_local(orb, 1, spin + 1, lat + 1) = Delta_local(orb, 1, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_1(ky))
          Delta_local(orb, 2, spin + 1, lat + 1) = Delta_local(orb, 2, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_2(kx, ky))
          Delta_local(orb, 3, spin + 1, lat + 1) = Delta_local(orb, 3, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_3(kx, ky))
        END DO
      END DO
    END DO
    !### END OF NEAREST NEIGHBOURS PAIRING ###########################################

    !### NEXT NEAREST NEIGHBOURS PAIRING ############################################
    !Electrons
    DO n = 1, DIM_POSITIVE_K
      occupation = fd_distribution(Energies(n), 0d0, T)
      !Up - down Ti1 - Ti1 delta, Ti2 - Ti2 delta
      !No conjugation in phase factor, since next nearest neighbours have the same relative positions in both sublattices
      DO lat = 0, SUBLATTICES - 1
        DO spin = 0, 1
          row = orb + lat * ORBITALS + DIM_POSITIVE_K + TBA_DIM * spin
          col = orb + lat * ORBITALS + TBA_DIM * MOD(spin + 1, 2)
          Delta_local(orb, N_NEIGHBOURS + 1, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 1, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_1(kx))
          Delta_local(orb, N_NEIGHBOURS + 2, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 2, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_2(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 3, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 3, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_3(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 4, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 4, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_4(kx))
          Delta_local(orb, N_NEIGHBOURS + 5, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 5, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_5(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 6, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 6, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_6(kx, ky))
        END DO
      END DO
    END DO

    !Holes
    DO n = DIM_POSITIVE_K + 1, DIM
      occupation = 1.-fd_distribution(-Energies(n), 0d0, T)
      !Up - down Ti1 - Ti1 delta, Ti2 - Ti2 delta
      !No conjugation in phase factor, since next nearest neighbours have the same relative positions in both sublattices
      DO lat = 0, SUBLATTICES - 1
        DO spin = 0, 1
          row = orb + lat * ORBITALS + DIM_POSITIVE_K + TBA_DIM * spin
          col = orb + lat * ORBITALS + TBA_DIM * MOD(spin + 1, 2)
          Delta_local(orb, N_NEIGHBOURS + 1, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 1, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_1(kx))
          Delta_local(orb, N_NEIGHBOURS + 2, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 2, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_2(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 3, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 3, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_3(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 4, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 4, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_4(kx))
          Delta_local(orb, N_NEIGHBOURS + 5, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 5, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_5(kx, ky))
          Delta_local(orb, N_NEIGHBOURS + 6, spin + 1, lat + 1) = Delta_local(orb, N_NEIGHBOURS + 6, spin + 1, lat + 1) + CONJG(U_transformation(row, n)) * U_transformation(col, n) * occupation * CONJG(pairing_nnn_6(kx, ky))
        END DO
      END DO
    END DO
    !### END OF NEXT NEAREST NEIGHBOURS PAIRING #####################################

  END DO

  !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
  Charge_dens_local(:) = 0.
  !Charge density calculation
  DO m = 1, DIM_POSITIVE_K
    DO n = 1, DIM_POSITIVE_K
      Charge_dens_local(m) = Charge_dens_local(m) + REAL(U_transformation(m, n) * CONJG(U_transformation(m, n)), KIND=8) * fd_distribution(Energies(n), 0d0, T) + &
      & REAL(U_transformation(m, DIM_POSITIVE_K + n) * CONJG(U_transformation(m, DIM_POSITIVE_K + n)), KIND=8) * (1.-fd_distribution(-Energies(DIM_POSITIVE_K + n), 0d0, T))
    END DO
  END DO

  !Multiplication by the Jacobian
  Delta_local = Delta_local * k1
  Charge_dens_local = Charge_dens_local * k1
END SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA

END MODULE mod_local_integrand
