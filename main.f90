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
    COMPLEX*16, ALLOCATABLE :: Delta(:,:,:,:), Delta_new(:,:,:,:)
    REAL*8, ALLOCATABLE :: Delta_broyden(:), Delta_new_broyden(:)
    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:), Gamma_SC_new(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    REAL*8 :: kx, ky

    INTEGER*4 :: i,j,n, i_orb, lat, a,b,c, m, orb, orb_prime,spin
    INTEGER*4 :: sc_iter
    LOGICAL :: sc_flag

    INTEGER*4 :: counter
    INTEGER*4 :: delta_real_elems
    INTEGER*4 :: broyden_index

    REAL*8 :: filling, filling_total

    !3 because of neighbours, 
    !2 because spin up-down and down-up,
    !2 because of complex number
    !SUBLATTICES because of Ti1-Ti2 coupling and Ti2 - Ti1 coupling (stored in this order)
    delta_real_elems = ORBITALS*3*2*2*SUBLATTICES

    CALL GET_INPUT("./input.nml")

    !Basis 
    !c_{k,yz,Ti1,up}, c_{k,zx,Ti1,up}, c_{k,xy,Ti1,up},
    !c_{k,yz,Ti2,up}, c_{k,zx,Ti2,up}, c_{k,xy,Ti2,up},
    !c_{k,yz,Ti1,down}, c_{k,zx,Ti1,down}, c_{k,xy,Ti1,down},
    !c_{k,yz,Ti2,down}, c_{k,zx,Ti2,down}, c_{k,xy,Ti2,down},
    ! + H.c{-k}
    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM,DIM))
    ALLOCATE(Energies(DIM))
    ALLOCATE(Delta(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))    !Third dimension for spin coupling up-down and down-up
    ALLOCATE(Delta_new(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Gamma_SC(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Gamma_SC_new(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Delta_broyden(delta_real_elems))   !Flattened Gamma array
    ALLOCATE(Delta_new_broyden(delta_real_elems))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))

    !Initializations
    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0., 0.)
    U_transformation(:,:) = DCMPLX(0., 0.)
    Energies(:) = 0.

    Delta(:,:,:,:) = DCMPLX(0. , 0.)
    Delta_new(:,:,:,:) = DCMPLX(0. , 0.)
    Gamma_SC(:,:,1,:) = DCMPLX(gamma_start, 0.)
    Gamma_SC(:,:,2,:) = DCMPLX(-gamma_start, 0.)
    !Gamma_SC(:,:,:) = DCMPLX(0., 0.)
    Gamma_SC_new(:,:,:,:) = DCMPLX(0., 0.)


    Delta_broyden(:) = 0.
    Delta_new_broyden(:) = 0.
    
    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:)) !This is not needed, since ZHEEV takes only upper triangle
    
    OPEN(unit = 99, FILE= "./OutputData/Convergence.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO sc_iter = 1, max_sc_iter
        !PRINT*, "============= SC_ITER: ", sc_iter
        !PRINT*, "Gamma = ", Gamma_SC(1,1,1,1)/meV2au
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
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1

                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))

                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:)) !This is not needed, since ZHEEV takes only upper triangle

                Hamiltonian(:,:) = 0.5*( Hamiltonian_const(:,:) + Hamiltonian(:,:) )
                !U_transformation(:,:) = Hamiltonian(:,:)
                !CALL DIAGONALIZE_HERMITIAN(U_transformation(:,:), Energies(i,j,:))
                !CALL PRINT_HAMILTONIAN(Hamiltonian(:,:))

                CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:))
                !After DIAGONALIZE HERMITIAN, U contains eigenvectors, so it corresponds to transformation matrix U                

                !Self - consistent delta calculation
                DO orb = 1, ORBITALS
                    DO lat = 0, SUBLATTICES - 1
                        !Electrons
                        DO n = 1, DIM_POSITIVE_K
                            !Up - down Ti1 - Ti2 delta
                            Delta_new(orb,1,1,1) = Delta_new(orb,1,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
                            Delta_new(orb,2,1,1) = Delta_new(orb,2,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx,ky)
                            Delta_new(orb,3,1,1) = Delta_new(orb,3,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx,ky)

                            !Up - down Ti2 - Ti1 delta
                            Delta_new(orb,1,1,2) = Delta_new(orb,1,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_1(ky))
                            Delta_new(orb,2,1,2) = Delta_new(orb,2,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_2(kx,ky))
                            Delta_new(orb,3,1,2) = Delta_new(orb,3,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_3(kx,ky))

                            !Down - up Ti1 - Ti2 delta
                            Delta_new(orb,1,2,1) = Delta_new(orb,1,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_1(ky)
                            Delta_new(orb,2,2,1) = Delta_new(orb,2,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_2(kx,ky)
                            Delta_new(orb,3,2,1) = Delta_new(orb,3,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*fd_distribution(Energies(n), 0d0, T)*pairing_3(kx,ky)
                            
                            !Down - up Ti2 - Ti1 delta
                            Delta_new(orb,1,2,2) = Delta_new(orb,1,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_1(ky))
                            Delta_new(orb,2,2,2) = Delta_new(orb,2,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_2(kx,ky))
                            Delta_new(orb,3,2,2) = Delta_new(orb,3,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*fd_distribution(Energies(n), 0d0, T)*CONJG(pairing_3(kx,ky))


                            !PRINT*, fd_distribution(Energies(i,j,n), 0d0, T)
                            ! filling = filling + fd_distribution(Energies(n), 0d0, T)
                            filling_total = filling_total + 1
                        END DO

                        !Holes
                        DO n = DIM_POSITIVE_K + 1, DIM
                            !Up - down Ti1 - Ti2 delta
                            Delta_new(orb,1,1,1) = Delta_new(orb,1,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
                            Delta_new(orb,2,1,1) = Delta_new(orb,2,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx,ky)
                            Delta_new(orb,3,1,1) = Delta_new(orb,3,1,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K, n))*U_transformation(orb + ORBITALS + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx,ky)

                            !Up - down Ti2 - Ti1 delta
                            Delta_new(orb,1,1,2) = Delta_new(orb,1,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_1(ky))
                            Delta_new(orb,2,1,2) = Delta_new(orb,2,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_2(kx,ky))
                            Delta_new(orb,3,1,2) = Delta_new(orb,3,1,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K, n))*U_transformation(orb + TBA_DIM, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_3(kx,ky))

                            !Down - up Ti1 - Ti2 delta
                            Delta_new(orb,1,2,1) = Delta_new(orb,1,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_1(ky)
                            Delta_new(orb,2,2,1) = Delta_new(orb,2,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_2(kx,ky)
                            Delta_new(orb,3,2,1) = Delta_new(orb,3,2,1) + CONJG(U_transformation(orb + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb + ORBITALS, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*pairing_3(kx,ky)
                            
                            !Down - up Ti2 - Ti1 delta
                            Delta_new(orb,1,2,2) = Delta_new(orb,1,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_1(ky))
                            Delta_new(orb,2,2,2) = Delta_new(orb,2,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_2(kx,ky))
                            Delta_new(orb,3,2,2) = Delta_new(orb,3,2,2) + CONJG(U_transformation(orb + ORBITALS + DIM_POSITIVE_K + TBA_DIM, n))*U_transformation(orb, n)*(1. - fd_distribution(-Energies(n), 0d0, T))*CONJG(pairing_3(kx,ky))

                            !PRINT*, 1. - fd_distribution(-Energies(i,j,n), 0d0, T)
                            ! filling = filling + (1. - fd_distribution(Energies(n), 0d0, T))
                            filling_total = filling_total + 1
                        END DO                        
                    END DO
                END DO

                DO m = 1, DIM_POSITIVE_K
                    DO n = 1, DIM_POSITIVE_K
                        Charge_dens(m) = Charge_dens(m) + ABS(U_transformation(m,n))**2 * fd_distribution(Energies(n), 0d0, T) + &
                        & ABS(U_transformation(m, DIM_POSITIVE_K + n))**2 * (1. - fd_distribution(-Energies(DIM_POSITIVE_K + n), 0d0, T))
                    END DO
                END DO


            END DO
        END DO !End of k-loop

        Delta_new(:,:,:,:) = Delta_new(:,:,:,:) * domega ! because of changing sum to integral
        Charge_dens(:) = Charge_dens(:) * domega
        !Gamma calculation
        DO spin = 1,2   !Loop over spin coupling up-down or down-up
            DO n = 1, N_NEIGHBOURS !Loop over neighbours
                DO lat = 1, SUBLATTICES !Loop over sublattices coupling
                    DO orb = 1, ORBITALS
                        Gamma_SC_new(orb,n,spin,lat) = -0.5*J_SC*Delta_new(orb,n,spin,lat)
                        DO orb_prime = 1, ORBITALS
                            IF(orb .NE. orb_prime) THEN
                                ! CHECK WHETHER THIS IS 0.25 OR 0.5
                                Gamma_SC_new(orb,n,spin,lat) = Gamma_SC_new(orb,n,spin,lat) - 0.5 * J_SC_PRIME * Delta_new(orb_prime, n,spin,lat)
                            END IF
                        END DO
                    END DO
                END DO
            END DO
        END DO

        !PRINT*, "Gamma new = ", Gamma_SC_new(1,1,1,1)/meV2au
        !PRINT*, "Filling ", SUM(Charge_dens(:))/(domega*k1_steps*k2_steps) / DIM_POSITIVE_K
        filling = 0
        filling_total = 0

        WRITE(99,'(I0, 4E15.5)') sc_iter, REAL(Gamma_SC(1,1,1,1)/meV2au), AIMAG(Gamma_SC(1,1,1,1)/meV2au), &
        &                                 REAL(Gamma_SC_new(1,1,1,1)/meV2au), AIMAG(Gamma_SC_new(1,1,1,1)/meV2au)

        
        !Here we should check whether convergence was reached
        sc_flag = .TRUE.
        DO spin = 1, 2
            DO orb = 1, ORBITALS
                DO n = 1, N_NEIGHBOURS
                    DO lat = 1, SUBLATTICES
                        ! IF ( ( ABS(REAL(Delta_new(orb,m,n) - Delta(orb,m,n))*J_SC) > eps_convergence ) .OR. &
                        ! & ( ABS(AIMAG(Delta_new(orb,m,n) - Delta(orb,m,n))*J_SC) > eps_convergence ) ) THEN
                        ! PRINT*, ABS( (ABS(Delta_new(orb,m,n)) - ABS(Delta(orb,m,n)) )*J_SC )/meV2au, eps_convergence/meV2au
                        IF(ABS( ABS(Gamma_SC_new(orb,n,spin,lat)) - ABS(Gamma_SC(orb,n,spin,lat)) )  > eps_convergence) THEN                   
                            sc_flag = .FALSE.
                        EXIT !Maybe go to???
                        END IF
                    END DO
                END DO
            END DO
        END DO

        IF (sc_flag) THEN 
            PRINT*, "Convergence reached!"
            EXIT
        END IF

        !Broyden mixing
        !Flatten arrays for Broyden mixing
        broyden_index = 1
        DO spin = 1, 2
            DO orb = 1, ORBITALS
                DO n = 1, N_NEIGHBOURS
                    DO lat = 1, SUBLATTICES
                        Delta_broyden(broyden_index) = REAL(Gamma_SC(orb,n,spin,lat))
                        Delta_broyden(INT(delta_real_elems/2) + broyden_index) = AIMAG(Gamma_SC(orb,n,spin,lat))
                        Delta_new_broyden(broyden_index) = REAL(Gamma_SC_new(orb,n,spin,lat))
                        Delta_new_broyden(INT(delta_real_elems/2) + broyden_index) = AIMAG(Gamma_SC_new(orb,n,spin,lat))
                        broyden_index = broyden_index + 1
                    END DO
                END DO
            END DO
        END DO

        CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .FALSE.)        

        broyden_index = 1
        DO spin = 1, 2
            DO orb = 1, ORBITALS
                DO n = 1, N_NEIGHBOURS
                    DO lat = 1, SUBLATTICES
                        Gamma_SC(orb,n,spin,lat) = DCMPLX(Delta_broyden(broyden_index), Delta_broyden(INT(delta_real_elems/2) + broyden_index))
                        broyden_index = broyden_index + 1
                    END DO
                END DO
            END DO
        END DO

        !Linear mixing
        !Gamma_SC(:,:,:) = (1. - sc_alpha)*Gamma_SC(:,:,:) + sc_alpha*Gamma_SC_new(:,:,:)

        Delta_new(:,:,:,:) = DCMPLX(0. , 0.)
        Gamma_SC_new(:,:,:,:) = DCMPLX(0., 0.)

    END DO !End of SC loop
    CLOSE(99)


    CALL PRINT_GAMMA(Gamma_SC(:,:,:,:), "Gamma_SC")




    !Just for memory deallocation
    CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .TRUE.)

    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(Energies)
    DEALLOCATE(Delta)
    DEALLOCATE(Delta_new)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Gamma_SC_new)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Delta_broyden)
    DEALLOCATE(Delta_new_broyden)
    DEALLOCATE(Charge_dens)

END PROGRAM MAIN