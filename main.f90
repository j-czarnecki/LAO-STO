PROGRAM MAIN
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_broydenV2
    USE mod_compute_hamiltonians
    USE mod_integrate

    IMPLICIT NONE 

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:)
    COMPLEX*16, ALLOCATABLE :: Delta_local(:,:,:,:), Delta_new(:,:,:,:)
    REAL*8, ALLOCATABLE :: Delta_broyden(:), Delta_new_broyden(:)
    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:), Gamma_SC_new(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:), Charge_dens_new(:), Charge_dens_local(:)

    REAL*8 :: gamma_error, gamma_max_error, charge_error, charge_max_error

    INTEGER*4 :: i,j,n, lat, orb, orb_prime,spin
    INTEGER*4 :: sc_iter
    LOGICAL :: sc_flag

    INTEGER*4 :: counter
    INTEGER*4 :: delta_real_elems
    INTEGER*4 :: broyden_index

    !3 because of neighbours, 
    !2 because spin up-down and down-up,
    !2 because of complex number
    !SUBLATTICES because of Ti1-Ti2 coupling and Ti2 - Ti1 coupling (stored in this order)
    !DIM_POSITIVE_K included due to Charge density self-consistency
    delta_real_elems = DIM_POSITIVE_K + ORBITALS*3*2*2!*SUBLATTICES !SUBLATTICES should be excluded in absence of magnetic field

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
    ALLOCATE(Delta_local(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))    !Third dimension for spin coupling up-down and down-up
    ALLOCATE(Delta_new(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Gamma_SC(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Gamma_SC_new(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Delta_broyden(delta_real_elems))   !Flattened Gamma array
    ALLOCATE(Delta_new_broyden(delta_real_elems))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))
    ALLOCATE(Charge_dens_local(DIM_POSITIVE_K))
    ALLOCATE(Charge_dens_new(DIM_POSITIVE_K))


    !Initializations
    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0., 0.)
    U_transformation(:,:) = DCMPLX(0., 0.)
    Energies(:) = 0.

    Delta_local(:,:,:,:) = DCMPLX(0. , 0.)
    Delta_new(:,:,:,:) = DCMPLX(0. , 0.)

    !Breaking spin up-down down-up symmetry
    Gamma_SC(:,:,1,:) = DCMPLX(gamma_start, 0.)
    Gamma_SC(:,:,2,:) = DCMPLX(-gamma_start, 0.)


    !Gamma_SC(:,:,:) = DCMPLX(0., 0.)
    Gamma_SC_new(:,:,:,:) = DCMPLX(0., 0.)

    Charge_dens(:) = charge_start
    Charge_dens_new(:) = 0.
    Charge_dens_local(:) = 0.

    Delta_broyden(:) = 0.
    Delta_new_broyden(:) = 0.
    
    gamma_max_error = 0.
    charge_max_error = 0.

    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
    
    OPEN(unit = 99, FILE= "./OutputData/Convergence.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO sc_iter = 1, max_sc_iter
        !PRINT*, "============= SC_ITER: ", sc_iter
        !PRINT*, "Gamma = ", Gamma_SC(1,1,1,1)/meV2au
        counter = 0

        !Those loops are only for slicing and future parallelization
        !Integration over chunks is computed via Romberg algorithm.
        DO i = 0, k1_steps-1
            DO j = 0, k2_steps-1    
                counter = counter + 1
                !PRINT*, counter
                
                CALL ROMBERG_Y(Hamiltonian_const(:,:), Gamma_SC(:,:,:,:), Charge_dens(:), i*dk1, (i + 1)*dk1, j*dk2, (j + 1)*dk2, &
                & Delta_local(:,:,:,:), Charge_dens_local(:), romb_eps_x, interpolation_deg_x, max_grid_refinements_x, &
                & romb_eps_y, interpolation_deg_y, max_grid_refinements_y)

                !This has to be atomic operations, since Delta_new and Charge_dens would be global variables for all threads
                Delta_new(:,:,:,:) = Delta_new(:,:,:,:) + Delta_local(:,:,:,:)
                Charge_dens_new(:) = Charge_dens_new(:) + Charge_dens_local(:)
            END DO
        END DO !End of k-loop
        !Due to change of sum to integral one has to divide by Brillouin zone volume
        Delta_new(:,:,:,:) = Delta_new(:,:,:,:)/(K1_MAX*K2_MAX)
        Charge_dens_new(:) = Charge_dens_new(:)/(K1_MAX*K2_MAX)

        ! !Integration with simple trapezoid rule
        ! DO i = 0, k1_steps
        !     counter = counter + 1
        !     PRINT*, counter
            
        !     DO j = 0, k2_steps

        !         CALL GET_LOCAL_CHARGE_AND_DELTA(Hamiltonian_const(:,:), Gamma_SC(:,:,:,:), &
        !         & Charge_dens(:), i*dk1, j*dk2, Delta_local(:,:,:,:), Charge_dens_local(:))
            
        !         Delta_new(:,:,:,:) = Delta_new(:,:,:,:) + Delta_local(:,:,:,:)*dk1*dk2
        !         Charge_dens_new(:) = Charge_dens_new(:) + Charge_dens_local(:)*dk1*dk2

        !     END DO
        ! END DO
        ! Delta_new(:,:,:,:) = Delta_new(:,:,:,:)/(K1_MAX*K2_MAX)
        ! Charge_dens_new(:) = Charge_dens_new(:)/(K1_MAX*K2_MAX)

        ! DO i = 1, DIM_POSITIVE_K
        !     PRINT*, "Charge elem ", i, " = ", Charge_dens_new(i)
        ! END DO


        !#########################################################################################################################
        !This is a critical section - only one thread can execute that and all thread should have ended their job up to that point
        !#########################################################################################################################

        !Delta_new(:,:,:,:) = Delta_new(:,:,:,:) * domega ! because of changing sum to integral
        !Charge_dens(:) = Charge_dens(:) * domega
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
        !PRINT*, "Filling ", SUM(Charge_dens(:)) / DIM_POSITIVE_K 
        
        !Here we check whether convergence was reached
        sc_flag = .TRUE.
        DO spin = 1, 2
            DO orb = 1, ORBITALS
                DO n = 1, N_NEIGHBOURS
                    DO lat = 1, 1 !to 1 in absence of magnetic field to SUBLATTICES if else
                        !It should be considered whether relative or absolute error must be checked
                        gamma_error = ABS( ABS(Gamma_SC_new(orb,n,spin,lat)) - ABS(Gamma_SC(orb,n,spin,lat)) )
                        !Gamma convergence checking
                        IF (gamma_error > gamma_eps_convergence) THEN                   
                            sc_flag = .FALSE.
                            !EXIT !Maybe go to???
                        END IF
                        
                        !Find biggest error in current iteration
                        IF (gamma_error > gamma_max_error) gamma_max_error = gamma_error

                    END DO
                END DO
            END DO
        END DO

        DO n = 1, DIM_POSITIVE_K
            charge_error = ABS(Charge_dens(n) - Charge_dens_new(n))
            IF (charge_error > charge_eps_convergence) THEN
                sc_flag = .FALSE.
            END IF

            IF (charge_error > charge_max_error) charge_max_error = charge_error
        END DO

        IF (sc_flag) THEN 
            PRINT*, "Convergence reached!"
            EXIT
        END IF

        WRITE(99,'(I0, 8E15.5)') sc_iter, REAL(Gamma_SC(1,1,1,1)/meV2au), AIMAG(Gamma_SC(1,1,1,1)/meV2au), &
        &                                 REAL(Gamma_SC_new(1,1,1,1)/meV2au), AIMAG(Gamma_SC_new(1,1,1,1)/meV2au), &
        &                                 Charge_dens(1), Charge_dens_new(1), gamma_max_error/meV2au, charge_max_error
        !PRINT*, "Gamma max error ", gamma_max_error

        !In the beginning of convergence use Broyden method to quickly find minimum
        ! IF (gamma_max_error  > 1e-3) THEN 
            !PRINT*, "Broyden mixing"
            !Broyden mixing
            !Flatten arrays for Broyden mixing
            broyden_index = 1
            DO spin = 1, 2
                DO orb = 1, ORBITALS
                    DO n = 1, N_NEIGHBOURS
                        DO lat = 1, 1 !to 1 in absence of magnetic field to SUBLATTICES if else
                            Delta_broyden(broyden_index) = REAL(Gamma_SC(orb,n,spin,lat))
                            Delta_broyden(INT((delta_real_elems - DIM_POSITIVE_K)/2) + broyden_index) = AIMAG(Gamma_SC(orb,n,spin,lat))
                            Delta_new_broyden(broyden_index) = REAL(Gamma_SC_new(orb,n,spin,lat))
                            Delta_new_broyden(INT((delta_real_elems - DIM_POSITIVE_K)/2) + broyden_index) = AIMAG(Gamma_SC_new(orb,n,spin,lat))
                            broyden_index = broyden_index + 1
                        END DO
                    END DO
                END DO
            END DO

            !Must be +1!!!
            Delta_broyden((delta_real_elems - DIM_POSITIVE_K + 1) : delta_real_elems) = Charge_dens(:)
            Delta_new_broyden((delta_real_elems - DIM_POSITIVE_K + 1) : delta_real_elems) = Charge_dens_new(:)

            !PRINT*, "Filled table for broyden mixing, calling mix_broyden"
            CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .FALSE.)        

            !PRINT*, "Finished mix_broyden, rewriting to Gamma"
            broyden_index = 1
            DO spin = 1, 2
                DO orb = 1, ORBITALS
                    DO n = 1, N_NEIGHBOURS
                        DO lat = 1, 1 !to 1 in absence of magnetic field to SUBLATTICES if else
                            Gamma_SC(orb,n,spin,lat) = DCMPLX(Delta_broyden(broyden_index), Delta_broyden(INT((delta_real_elems - DIM_POSITIVE_K)/2) + broyden_index))
                            broyden_index = broyden_index + 1
                        END DO
                    END DO
                END DO
            END DO

            Charge_dens(:) = Delta_broyden((delta_real_elems - DIM_POSITIVE_K + 1):delta_real_elems)
        !In the last phase of convergence use linear mixing to avoid spare oscillations
        ! ELSE
        !     PRINT*, "Linear mixing"
        !     !Linear mixing
        !     Gamma_SC(:,:,:,:) = (1. - sc_alpha)*Gamma_SC(:,:,:,:) + sc_alpha*Gamma_SC_new(:,:,:,:)
        ! END IF


        Gamma_SC(:,:,:,2) = CONJG(Gamma_SC(:,:,:,1)) !This is valid in asbence of magnetic field
        Delta_new(:,:,:,:) = DCMPLX(0. , 0.)
        Gamma_SC_new(:,:,:,:) = DCMPLX(0., 0.)
        Charge_dens_new(:) = 0.
        gamma_max_error = 0.
        charge_max_error = 0.

        !To check the state of the simulation
        CALL PRINT_GAMMA(Gamma_SC(:,:,:,:), "Gamma_SC_iter")
        CALL PRINT_CHARGE(Charge_dens(:), "Chargen_dens_iter")

    END DO !End of SC loop
    CLOSE(99)

    !Printing results after the simulation is done
    CALL PRINT_GAMMA(Gamma_SC(:,:,:,:), "Gamma_SC_final")
    CALL PRINT_CHARGE(Charge_dens(:), "Chargen_dens_final")


    !Just for memory deallocation, the .TRUE. flag is crucial
    CALL mix_broyden(delta_real_elems, Delta_new_broyden(:), Delta_broyden(:), sc_alpha, sc_iter, 4, .TRUE.)




    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(Energies)
    DEALLOCATE(Delta_local)
    DEALLOCATE(Delta_new)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Gamma_SC_new)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Delta_broyden)
    DEALLOCATE(Delta_new_broyden)
    DEALLOCATE(Charge_dens)
    DEALLOCATE(Charge_dens_local)
    DEALLOCATE(Charge_dens_new)

END PROGRAM MAIN