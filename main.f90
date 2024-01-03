PROGRAM MAIN
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_broydenV2

    IMPLICIT NONE 

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:)
    COMPLEX*16, ALLOCATABLE :: Delta(:,:,:), Delta_new(:,:,:)
    REAL*8, ALLOCATABLE :: Delta_broyden(:), Delta_new_broyden(:)

    REAL*8 :: kx, ky

    INTEGER*4 :: i,j,n, i_orb, lat, a,b,c, m, orb
    INTEGER*4 :: sc_iter
    LOGICAL :: sc_flag

    INTEGER*4 :: counter
    INTEGER*4 :: delta_real_elems
    INTEGER*4 :: broyden_index

    REAL*8 :: filling, filling_total

    !3 because of neghbours, 
    !2 because spin up-down and down-up,
    !2 because of complex number 
    delta_real_elems = ORBITALS*3*2*2 

    CALL GET_INPUT("./input.nml")

    ALLOCATE(Hamiltonian(DIM,DIM))
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM,DIM))
    ALLOCATE(Energies(DIM))
    ALLOCATE(Delta(ORBITALS, N_NEIGHBOURS,2))   !Third dimension for spin coupling: 1 - up-down coupling, 2 - down-up coupling
    ALLOCATE(Delta_new(ORBITALS, N_NEIGHBOURS,2))
    ALLOCATE(Delta_broyden(delta_real_elems))
    ALLOCATE(Delta_new_broyden(delta_real_elems))

    !Initializations
    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0., 0.)
    U_transformation(:,:) = DCMPLX(0., 0.)
    Energies(:) = 0.

    Delta(:,:,:) = DCMPLX(1e-3 , 1e-3)
    !Delta(:,:,:) = DCMPLX(0. , 0.)
    Delta_new(:,:,:) = DCMPLX(0. , 0.)

    Delta_broyden(:) = 0.
    Delta_new_broyden(:) = 0.
    
    !Computing k-independent terms
    ! CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    ! CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    ! CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:)) !This is not needed, since ZHEEV takes only upper triangle
    
    OPEN(unit = 99, FILE= "./OutputData/Convergence.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO sc_iter = 1, max_sc_iter
        !PRINT*, "============= SC_ITER: ", sc_iter
        OPEN(unit = 9, FILE= "./OutputData/Ek.dat", FORM = "FORMATTED", ACTION = "WRITE")

        counter = 0
        DO i = 0, k1_steps
        !DO i = -k1_steps,k1_steps
            !DO j = -ky_steps, ky_steps
            !DO j = -MIN(k1_steps - i, k1_steps), MIN(k1_steps + i, k1_steps) !This guarantees integrating over first Brillouin zone
            DO j = 0, k2_steps    
                counter = counter + 1
                ! kx = i*dkx*A_TILDE
                ! ky = j*dky*A_TILDE
                !Transform from graphene reciprocal lattice to kx and ky
                ! kx = i*dk1/SQRT(3.) * A_TILDE
                ! ky = (-i*dk1/3. + 2./3.*j*dk2) * A_TILDE
                kx = ( i*dk1*SQRT(3.)/2. ) * A_TILDE
                ky = ( -i*dk1/2. + j*dk2 ) * A_TILDE
                Energies(:) = 0.
                Hamiltonian(:,:) = DCMPLX(0. , 0.)
                U_transformation(:,:) = DCMPLX(0. , 0.)
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                ! CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)
                ! CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky)
                ! CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)

                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Delta(:,:,:))

                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:)) !This is not needed, since ZHEEV takes only upper triangle

                Hamiltonian(:,:) = 0.5*( Hamiltonian_const(:,:) + Hamiltonian(:,:) )
                !U_transformation(:,:) = Hamiltonian(:,:)
                !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:))
                CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:))
                !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U                
                DO n = 1, DIM
                    WRITE(9, *) kx, ky, Energies(n)
                END DO

                !Self - consistent delta calculation
                DO orb = 1, ORBITALS
                    DO lat = 0, SUBLATTICES - 1

                        !Energies are sorted from lowest to highest
                        !Electrons
                        DO n = 1, DIM_POSITIVE_K
                            !Up - down coupling
                            Delta_new(orb,1,1) = Delta_new(orb,1,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
                            Delta_new(orb,2,1) = Delta_new(orb,2,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx, ky)
                            Delta_new(orb,3,1) = Delta_new(orb,3,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx, ky)

                            !Down - up coupling
                            Delta_new(orb,1,2) = Delta_new(orb,1,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
                            Delta_new(orb,2,2) = Delta_new(orb,2,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx, ky)
                            Delta_new(orb,3,2) = Delta_new(orb,3,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx, ky)
                            !PRINT*, fd_distribution(Energies(i,j,n), 0d0, T)
                            filling = filling + fd_distribution(Energies(n), 0d0, T)
                            filling_total = filling_total + 1
                        END DO

                        !Holes
                        DO n = DIM_POSITIVE_K + 1, DIM
                            !Up - down coupling
                            Delta_new(orb,1,1) = Delta_new(orb,1,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
                            Delta_new(orb,2,1) = Delta_new(orb,2,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx, ky)
                            Delta_new(orb,3,1) = Delta_new(orb,3,1) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + lat*ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx, ky)

                            !Down - up coupling
                            Delta_new(orb,1,2) = Delta_new(orb,1,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
                            Delta_new(orb,2,2) = Delta_new(orb,2,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx, ky)
                            Delta_new(orb,3,2) = Delta_new(orb,3,2) + CONJG(U_transformation(orb + lat*ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + lat*ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx, ky)
                            !PRINT*, 1. - fd_distribution(-Energies(i,j,n), 0d0, T)
                            filling = filling + (1. - fd_distribution(Energies(n), 0d0, T))
                            filling_total = filling_total + 1
                        END DO                        
                    END DO
                END DO


            END DO
            WRITE(9,*)
            WRITE(9,*)
        END DO !End of k-loop
        CLOSE(9)


        !PRINT*, "Filling/filling_total ", filling/filling_total
        filling = 0
        filling_total = 0
        !PRINT*, "N sites ", counter
        Delta_new = Delta_new * domega ! * V/(2*PI)**2, because of changing sum to integral?
        !Delta_new(:,:,1) = Delta_new(:,:,2) !For test: assuming that up-down coupling is the same as down-up
        ! PRINT*, "Gamma 1,1", Delta(1,1,1)*J_SC/meV2au
        ! PRINT*, "Abs Gamma 1,1", ABS(Delta(1,1,1))*J_SC/meV2au
        ! PRINT*, "Gamma new 1,1", Delta_new(1,1,1)*J_SC/meV2au
        ! PRINT*, "Abs Gamma new 1,1", ABS(Delta_new(1,1,1))*J_SC/meV2au

        WRITE(99,'(I0, 4E15.5)') sc_iter, REAL(Delta(1,1,1)*J_SC/meV2au), AIMAG(Delta(1,1,1)*J_SC/meV2au), &
        &                                 REAL(Delta_new(1,1,1)*J_SC/meV2au), AIMAG(Delta_new(1,1,1)*J_SC/meV2au)

        
        !Here we should check whether convergence was reached
        sc_flag = .TRUE.
        DO orb = 1, ORBITALS
            DO m = 1, 3
                DO n = 1, 2
                    ! IF ( ( ABS(REAL(Delta_new(orb,m,n) - Delta(orb,m,n))*J_SC) > eps_convergence ) .OR. &
                    ! & ( ABS(AIMAG(Delta_new(orb,m,n) - Delta(orb,m,n))*J_SC) > eps_convergence ) ) THEN
                    ! PRINT*, ABS( (ABS(Delta_new(orb,m,n)) - ABS(Delta(orb,m,n)) )*J_SC )/meV2au, eps_convergence/meV2au
                    IF(ABS( ( ABS(Delta_new(orb,m,n)) - ABS(Delta(orb,m,n)) )*J_SC )  > eps_convergence) THEN                   
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

        !Broyden mixing
        !Flatten arrays for Broyden mixing
        ! broyden_index = 1
        ! DO orb = 1, ORBITALS
        !     DO m = 1, 3
        !         DO n = 1, 2
        !             Delta_broyden(broyden_index) = REAL(Delta(orb,m,n))
        !             Delta_broyden(INT(delta_real_elems/2) + broyden_index) = AIMAG(Delta(orb,m,n))
        !             Delta_new_broyden(broyden_index) = REAL(Delta_new(orb,m,n))
        !             Delta_new_broyden(INT(delta_real_elems/2) + broyden_index) = AIMAG(Delta_new(orb,m,n))
        !             broyden_index = broyden_index + 1
        !         END DO
        !     END DO
        ! END DO

        ! CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 1, .FALSE.)        
        
        ! broyden_index = 1
        ! DO orb = 1, ORBITALS
        !     DO m = 1, 3
        !         DO n = 1, 2
        !             Delta(orb,m,n) = DCMPLX(Delta_broyden(broyden_index), Delta_broyden(INT(delta_real_elems/2) + broyden_index))
        !             broyden_index = broyden_index + 1
        !         END DO
        !     END DO
        ! END DO

        !Linear mixing
        Delta(:,:,:) = (1. - sc_alpha)*Delta(:,:,:) + sc_alpha*Delta_new(:,:,:)

        Delta_new(:,:,:) = DCMPLX(0. , 0.)


    END DO !End of SC loop
    CLOSE(99)


    CALL PRINT_DELTA(Delta(:,:,:), "Delta_SC")

    !Just for memory deallocation
    !CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 1, .TRUE.)

    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(Energies)
    DEALLOCATE(Delta)
    DEALLOCATE(Delta_new)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Delta_broyden)
    DEALLOCATE(Delta_new_broyden)

END PROGRAM MAIN