PROGRAM MAIN
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities

    IMPLICIT NONE 

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)
    COMPLEX*16, ALLOCATABLE :: Delta(:,:)

    REAL*8 :: kx, ky 

    REAL*8 :: dkx, dky, kx_max, ky_max !Suppose kx_max == ky_max and dk_x == dk_y
    INTEGER*4 :: kx_steps, ky_steps

    INTEGER*4 :: i,j,n
    INTEGER*4 :: sc_iter, max_sc_iter

    max_sc_iter = 1

    kx_max = 2. * 1./A_TILDE !Multiplier is arbitrary. Part of Brillouin zone to integrate over
    ky_max = 2. * 1./A_TILDE
    kx_steps = 200
    ky_steps = 200
    dkx = kx_max / kx_steps
    dky = ky_max / ky_steps

    ALLOCATE(Hamiltonian(DIM,DIM))
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(Energies(-kx_steps:kx_steps, -ky_steps:ky_steps, DIM))
    ALLOCATE(Delta(ORBITALS, N_NEIGHBOURS))

    Delta(:,:) = 1e3 * meV2au

    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:))

    DO sc_iter = 1, max_sc_iter
        IF (sc_iter .NE. 1) THEN
            !Here self-consistent superconducting deltas should be calculated
        END IF


        !OPEN(unit = 9, FILE= "./OutputData/E_k.dat", FORM = "FORMATTED", ACTION = "WRITE")
        DO i = -kx_steps,kx_steps
            DO j = -ky_steps, ky_steps
                kx = i*dkx*A_TILDE
                ky = j*dky*A_TILDE
                Energies(i,j,:) = 0.
                Hamiltonian(:,:) = 0.
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)

                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, delta(:,:))

                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:))

                Hamiltonian(:,:) = Hamiltonian_const(:,:) + Hamiltonian(:,:)
                Hamiltonian(:,:) = 1./2. * Hamiltonian(:,:)
                CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:,:), Energies(i,j,:))
                
                ! DO n = 1, DIM 
                !     WRITE(9,*) kx*A_TILDE, ky*A_TILDE, au_to_meV(Energies(i,j,n))*1e-3
                ! END DO 
                ! WRITE(9,*)
                ! WRITE(9,*)
            END DO
        END DO
        !CLOSE(9)


        !Here we should check whether convergence was reached

    END DO !End of SC loop

    OPEN(unit = 9, FILE= "./OutputData/E_k_sorted.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 1, DIM
        DO i = -kx_steps,kx_steps
            DO j = -ky_steps, ky_steps
                WRITE(9,*) i*dkx*A_TILDE, j*dky*A_TILDE, au_to_meV(Energies(i, j, n))*1e-3
            END DO
        END DO
        WRITE(9,*)
        WRITE(9,*)
    END DO 
    CLOSE(9)

    OPEN(unit = 9, FILE= "./OutputData/E_kx0_slice.dat", FORM = "FORMATTED", ACTION = "WRITE")
    OPEN(unit = 10, FILE= "./OutputData/E_ky0_slice.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 1, DIM
        DO i = 0, kx_steps
            WRITE(10, *) i*dkx*A_TILDE, au_to_meV(Energies(i,0,n))*1e-3
        END DO 
        DO i = -ky_steps, 0
            WRITE(9,*) i*dky*A_TILDE, au_to_meV(Energies(0, i, n))*1e-3
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

END PROGRAM MAIN