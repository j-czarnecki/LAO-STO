PROGRAM MAIN
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader

    IMPLICIT NONE 

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)
    COMPLEX*16, ALLOCATABLE :: Delta(:,:,:), Delta_new(:,:,:)

    REAL*8 :: kx, ky!, domega

    !REAL*8 :: dk1, dk2!, k1_max, k2_max !Suppose kx_max == ky_max and dk_x == dk_y
    !INTEGER*4 :: k1_steps, k2_steps

    INTEGER*4 :: i,j,n, i_orb, lat, a,b,c, m, orb
    INTEGER*4 :: sc_iter!, max_sc_iter
    !REAL*8 :: sc_alpha
    !REAL*8 :: eps_convergence
    LOGICAL :: sc_flag

    INTEGER*4 :: counter

    CALL GET_INPUT("./input.nml")

    !max_sc_iter = 1

    !k1_max = (2 * PI * 2./3.)/A_TILDE !Full Brillouin zone to integrate over
    !k2_max = (2 * PI * 2./3.)/A_TILDE
    ! kx_max = 2./A_TILDE
    ! ky_max = 2./A_TILDE
    !k1_steps = 500
    !k2_steps = 500
    !dk1 = k1_max / k1_steps
    !dk2 = k2_max / k2_steps
    !domega = ABS(dk1*dk2*SIN(PI/3.))

    !sc_alpha = 0.2
    !eps_convergence = 1e-5

    ALLOCATE(Hamiltonian(DIM,DIM))
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM,DIM))
    ALLOCATE(Energies(-k1_steps:k1_steps, -k2_steps:k2_steps, DIM))
    ALLOCATE(Delta(ORBITALS, N_NEIGHBOURS,2))   !Third dimension for spin coupling: 1 - up-down coupling, 2 - down-up coupling
    ALLOCATE(Delta_new(ORBITALS, N_NEIGHBOURS,2))

    Delta(:,:,:) = 0.
    Delta_new(:,:,:) = 0.
    
    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(n,n) + E_Fermi
    END DO
    !CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:)) !This is not needed, since ZHEEV takes only upper triangle

    DO sc_iter = 1, max_sc_iter
        PRINT*, "============= SC_ITER: ", sc_iter
        PRINT*, "Delta 1,1", Delta(1,1,1)

        counter = 0
        DO i = -k1_steps,k1_steps
            !DO j = -ky_steps, ky_steps
            DO j = -MIN(k1_steps - i, k1_steps), MIN(k1_steps + i, k1_steps) !This guarantees integrating over first Brillouin zone
                counter = counter + 1
                ! kx = i*dkx*A_TILDE
                ! ky = j*dky*A_TILDE
                !Transform from graphene reciprocal lattice to kx and ky
                kx = i*dk1/SQRT(3.) * A_TILDE
                ky = (-i*dk1/3. + 2./3.*j*dk2) * A_TILDE
                Energies(i,j,:) = 0.
                Hamiltonian(:,:) = 0.
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)

                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Delta(:,:,:))

                !CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:)) !This is not needed, since ZHEEV takes only upper triangle

                Hamiltonian(:,:) = 1./2.*( Hamiltonian_const(:,:) + Hamiltonian(:,:) )
                U_transformation(:,:) = Hamiltonian(:,:)
                !CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:,:), Energies(i,j,:))
                CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:))
                !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U                
                
                !CALL PRINT_HAMILTONIAN(REAL(Hamiltonian(:,:)), "H_real")
                !CALL PRINT_HAMILTONIAN(AIMAG(Hamiltonian(:,:)), "H_imag")
                ! CALL PRINT_HAMILTONIAN(REAL(U_transformation(:,:)), "U_transformation_real")
                ! CALL PRINT_HAMILTONIAN(AIMAG(U_transformation(:,:)), "U_transformation_imag")
                ! CALL PRINT_HAMILTONIAN(ABS(U_transformation(:,:))**2, "U_transformation_module")

                !Self - consistent delta calculation
                DO orb = 1, ORBITALS
                    DO lat = 0, SUBLATTICES - 1

                        !Energies are sorted from lowest to highest
                        !Holes
                        DO n = 1, DIM_POSITIVE_K
                            ! Up- down pairing
                            Delta_new(orb,1,1) = Delta_new(orb,1,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*( 1. - fd_distribution(Energies(i,j,n), 0d0 , T))*pairing_1(-ky) )
                            Delta_new(orb,2,1) = Delta_new(orb,2,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*( 1. - fd_distribution(Energies(i,j,n), 0d0 , T))*pairing_2(-kx, -ky) )
                            Delta_new(orb,3,1) = Delta_new(orb,3,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*( 1. - fd_distribution(Energies(i,j,n), 0d0 , T))*pairing_3(-kx, -ky) )
                        END DO

                        !Electrons
                        DO n = DIM_POSITIVE_K + 1, DIM
                            ! Up- down pairing
                            Delta_new(orb,1,1) = Delta_new(orb,1,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*fd_distribution(Energies(i,j,n), 0d0 , T)*pairing_1(-ky) )
                            Delta_new(orb,2,1) = Delta_new(orb,2,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*fd_distribution(Energies(i,j,n), 0d0 , T)*pairing_2(-kx, -ky) )
                            Delta_new(orb,3,1) = Delta_new(orb,3,1) + &
                            & CONJG( CONJG(U_transformation(n, orb + lat*SUBLATTICES))*U_transformation(orb + lat*SUBLATTICES + TBA_DIM + DIM_POSITIVE_K, n)*fd_distribution(Energies(i,j,n), 0d0 , T)*pairing_3(-kx, -ky) )
                        END DO                        
                    END DO
                END DO

            END DO
        END DO !End of k-loop

        PRINT*, "N sites ", counter
        Delta_new = Delta_new * domega
        Delta_new(:,:,2) = Delta_new(:,:,1) !For test: assuming that up-down coupling is the same as down-up
        !Here we should check whether convergence was reached
        PRINT*, "New delta 1,1", Delta_new(1,1,1)
        sc_flag = .TRUE.
        DO orb = 1, ORBITALS
            DO m = 1, 3
                DO n = 1, 2
                    IF ( (REAL(Delta_new(orb,m,n)) - REAL(Delta(orb,m,n)) > eps_convergence) .OR. &
                    & (AIMAG(Delta_new(orb,m,n)) - AIMAG(Delta(orb,m,n)) > eps_convergence)) THEN
                        sc_flag = .FALSE.
                        EXIT !Maybe go to???
                    END IF
                END DO
            END DO
        END DO

        IF (sc_flag) THEN 
            PRINT*, "Convergence reached!"
            EXIT
        END IF

        !Here self-consistent superconducting deltas should be calculated
        Delta(:,:,:) = (1. - sc_alpha)*Delta(:,:,:) + sc_alpha*Delta_new(:,:,:)
        Delta_new(:,:,:) = 0.


    END DO !End of SC loop

    CALL PRINT_ENERGIES(Energies(:,:,:), k1_steps, k2_steps, dk1, dk2, "Ek_sorted")
    CALL PRINT_DELTA(Delta(:,:,:), "Delta_SC")
    OPEN(unit = 9, FILE= "./OutputData/E_kx0_slice.dat", FORM = "FORMATTED", ACTION = "WRITE")
    OPEN(unit = 10, FILE= "./OutputData/E_ky0_slice.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 1, DIM
        DO i = 0, k1_steps
            WRITE(10, *) i*dk1 / k1_max, au_to_meV(Energies(i,0,n))*1e-3
        END DO 
        DO i = -k2_steps, 0
            WRITE(9,*) i*dk2 / k2_max, au_to_meV(Energies(0, i, n))*1e-3
        END DO

        WRITE(9,*)
        WRITE(9,*)
        WRITE(10,*)
        WRITE(10,*)
    END DO
    CLOSE(9)
    CLOSE(10)


    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(Energies)
    DEALLOCATE(Delta)
    DEALLOCATE(Delta_new)
    DEALLOCATE(U_transformation)

END PROGRAM MAIN