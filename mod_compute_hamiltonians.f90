MODULE mod_compute_hamiltonians
USE mod_parameters
USE mod_utilities
USE mod_hamiltonians
IMPLICIT NONE
CONTAINS

SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const, Gamma_SC, k1, k2, Delta_local, Charge_dens_local)
    COMPLEX*16, INTENT(IN) :: Hamiltonian_const(DIM, DIM)
    REAL*8, INTENT(IN) :: k1, k2
    COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES)

    COMPLEX*16, INTENT(OUT) :: Delta_local(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES)
    REAL*8, INTENT(OUT) :: Charge_dens_local(DIM_POSITIVE_K)
    
    COMPLEX*16 :: Hamiltonian(DIM, DIM)
    COMPLEX*16 :: U_transformation(DIM, DIM)
    REAL*8 :: Energies(DIM)
    REAL*8 :: kx, ky
    INTEGER*4 :: orb, lat, n, m


    !Transform from graphene reciprocal lattice to kx and ky
    kx = ( k1*SQRT(3.)/2. ) * A_TILDE
    ky = ( -k1/2. + k2 ) * A_TILDE
    Energies(:) = 0.
    Hamiltonian(:,:) = DCMPLX(0. , 0.)
    U_transformation(:,:) = DCMPLX(0. , 0.)
    CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
    CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
    CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens_local(:))
    CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))

    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:)) !This is not needed, since ZHEEV takes only upper triangle

    Hamiltonian(:,:) = 0.5*( Hamiltonian_const(:,:) + Hamiltonian(:,:) )
    !U_transformation(:,:) = Hamiltonian(:,:)
    !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:))
    !CALL PRINT_HAMILTONIAN(Hamiltonian(:,:))

    CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:))
    !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U                

    !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
    Delta_local(:,:,:,:) = DCMPLX(0. , 0.)
    !Self - consistent delta calculation
    DO orb = 1, ORBITALS
        DO lat = 0, SUBLATTICES - 1
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
        END DO
    END DO


    !Here it has to be set to zero, to avoid artifacts from previous iteration / chunk
    Charge_dens_local(:) = 0.
    !Charge density calculation
    DO m = 1, DIM_POSITIVE_K
        DO n = 1, DIM_POSITIVE_K
            Charge_dens_local(m) = Charge_dens_local(m) + ABS(U_transformation(m,n))**2 * fd_distribution(Energies(n), 0d0, T) + &
            & ABS(U_transformation(m, DIM_POSITIVE_K + n))**2 * (1. - fd_distribution(-Energies(DIM_POSITIVE_K + n), 0d0, T))
        END DO
    END DO

END SUBROUTINE GET_LOCAL_CHARGE_AND_DELTA



END MODULE mod_compute_hamiltonians