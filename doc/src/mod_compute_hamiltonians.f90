MODULE mod_compute_hamiltonians
USE mod_parameters
USE mod_utilities
USE mod_hamiltonians
USE mod_writers
IMPLICIT NONE
CONTAINS

SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, Charge_dens, k1, k2, Delta_local, Charge_dens_local)
    COMPLEX*16, INTENT(IN) :: Hamiltonian_const(DIM, DIM)
    REAL*8, INTENT(IN) :: k1, k2
    COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES)
    REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K)

    COMPLEX*16, INTENT(OUT) :: Delta_local(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES)
    REAL*8, INTENT(OUT) :: Charge_dens_local(DIM_POSITIVE_K)
    
    COMPLEX*16 :: Hamiltonian(DIM, DIM)
    COMPLEX*16 :: U_transformation(DIM, DIM)
    REAL*8 :: Energies(DIM)
    REAL*8 :: kx, ky
    INTEGER*4 :: orb, lat, n, m, spin


    !Transform from graphene reciprocal lattice to kx and ky
    kx = 2.*PI/(SQRT(3.0d0)) * k1
    ky = -2.*PI/3. * k1 + 4.*PI/3. * k2
    Energies(:) = 0.
    Hamiltonian(:,:) = DCMPLX(0. , 0.)
    U_transformation(:,:) = DCMPLX(0. , 0.)
    CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
    CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
    CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))

    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle

    Hamiltonian(:,:) = 0.5*( Hamiltonian_const(:,:) + Hamiltonian(:,:) )
    !U_transformation(:,:) = Hamiltonian(:,:)
    !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:), DIM)
    ! CALL PRINT_HAMILTONIAN(Hamiltonian(:,:))
    ! STOP 'Hamiltonian printed'

    CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), DIM)
    !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U                

    !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
    Delta_local(:,:,:,:) = DCMPLX(0. , 0.)
    !Self - consistent delta calculation
    DO orb = 1, ORBITALS
        !### NEAREST NEIGHBOURS PAIRING ###########################################
        !Electrons
        DO n = 1, DIM_POSITIVE_K
            !Up - down Ti1 - Ti2 delta
            Delta_local(orb,1,1,1) = Delta_local(orb,1,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
            Delta_local(orb,2,1,1) = Delta_local(orb,2,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx,ky)
            Delta_local(orb,3,1,1) = Delta_local(orb,3,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx,ky)

            !Up - down Ti2 - Ti1 delta
            Delta_local(orb,1,1,2) = Delta_local(orb,1,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_1(ky))
            Delta_local(orb,2,1,2) = Delta_local(orb,2,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_2(kx,ky))
            Delta_local(orb,3,1,2) = Delta_local(orb,3,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_3(kx,ky))

            !Down - up Ti1 - Ti2 delta
            Delta_local(orb,1,2,1) = Delta_local(orb,1,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
            Delta_local(orb,2,2,1) = Delta_local(orb,2,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx,ky)
            Delta_local(orb,3,2,1) = Delta_local(orb,3,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx,ky)
            
            !Down - up Ti2 - Ti1 delta
            Delta_local(orb,1,2,2) = Delta_local(orb,1,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_1(ky))
            Delta_local(orb,2,2,2) = Delta_local(orb,2,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_2(kx,ky))
            Delta_local(orb,3,2,2) = Delta_local(orb,3,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_3(kx,ky))
        END DO

        !Holes
        DO n = DIM_POSITIVE_K + 1, DIM
            !Up - down Ti1 - Ti2 delta
            Delta_local(orb,1,1,1) = Delta_local(orb,1,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
            Delta_local(orb,2,1,1) = Delta_local(orb,2,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx,ky)
            Delta_local(orb,3,1,1) = Delta_local(orb,3,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx,ky)

            !Up - down Ti2 - Ti1 delta
            Delta_local(orb,1,1,2) = Delta_local(orb,1,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_1(ky))
            Delta_local(orb,2,1,2) = Delta_local(orb,2,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_2(kx,ky))
            Delta_local(orb,3,1,2) = Delta_local(orb,3,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_3(kx,ky))

            !Down - up Ti1 - Ti2 delta
            Delta_local(orb,1,2,1) = Delta_local(orb,1,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
            Delta_local(orb,2,2,1) = Delta_local(orb,2,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx,ky)
            Delta_local(orb,3,2,1) = Delta_local(orb,3,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx,ky)
            
            !Down - up Ti2 - Ti1 delta
            Delta_local(orb,1,2,2) = Delta_local(orb,1,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_1(ky))
            Delta_local(orb,2,2,2) = Delta_local(orb,2,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_2(kx,ky))
            Delta_local(orb,3,2,2) = Delta_local(orb,3,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_3(kx,ky))

        END DO  
        !### END OF NEAREST NEIGHBOURS PAIRING ###########################################
                      
        !### NEXT NEAREST NEIGHBOURS PAIRING ############################################
        !Electrons
        DO n = 1, DIM_POSITIVE_K
            !Up - down Ti1 - Ti1 delta, Ti2 - Ti2 delta
            !No conjugation in phase factor, since next nearest neighbours have the same relative positions in both sublattices
            DO lat = 0, SUBLATTICES - 1
                DO spin = 0, 1
                    Delta_local(orb,N_NEIGHBOURS + 1,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 1,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_1(kx)
                    Delta_local(orb,N_NEIGHBOURS + 2,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 2,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_2(kx,ky)
                    Delta_local(orb,N_NEIGHBOURS + 3,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 3,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_3(kx, ky)
                    Delta_local(orb,N_NEIGHBOURS + 4,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 4,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_4(kx)
                    Delta_local(orb,N_NEIGHBOURS + 5,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 5,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_5(kx, ky)
                    Delta_local(orb,N_NEIGHBOURS + 6,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 6,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*fd_distribution(Energies(n), 0d0, T)*pairing_nnn_6(kx, ky)
                END DO
            END DO
        END DO

        !Holes
        DO n = DIM_POSITIVE_K + 1, DIM
            !Up - down Ti1 - Ti1 delta, Ti2 - Ti2 delta
            !No conjugation in phase factor, since next nearest neighbours have the same relative positions in both sublattices
            DO lat = 0, SUBLATTICES - 1
                DO spin = 0, 1
                    Delta_local(orb,N_NEIGHBOURS + 1,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 1,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_1(kx)
                    Delta_local(orb,N_NEIGHBOURS + 2,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 2,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_2(kx,ky)
                    Delta_local(orb,N_NEIGHBOURS + 3,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 3,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_3(kx, ky)
                    Delta_local(orb,N_NEIGHBOURS + 4,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 4,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_4(kx)
                    Delta_local(orb,N_NEIGHBOURS + 5,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 5,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_5(kx, ky)
                    Delta_local(orb,N_NEIGHBOURS + 6,spin + 1,lat + 1) = Delta_local(orb,N_NEIGHBOURS + 6,spin + 1,lat + 1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM*spin, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM*MOD(spin + 1, 2), n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_nnn_6(kx, ky)
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
            Charge_dens_local(m) = Charge_dens_local(m) + U_transformation(m,n)*CONJG(U_transformation(m,n)) * fd_distribution(Energies(n), 0d0, T) + &
            & U_transformation(m, DIM_POSITIVE_K + n)*CONJG(U_transformation(m, DIM_POSITIVE_K + n)) * (1. - fd_distribution(-Energies(DIM_POSITIVE_K + n), 0d0, T))
        END DO
    END DO

END SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA



END MODULE mod_compute_hamiltonians