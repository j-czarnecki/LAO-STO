PROGRAM dispersion
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_compute_hamiltonians
    IMPLICIT NONE

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    REAL*8, ALLOCATABLE :: DOS(:)

    REAL*8 :: k1, k2, kx, ky
    INTEGER*4 :: i,j,k,n, lat, orb, orb_prime,spin

    !DOS calculation
    REAL*8 :: E_DOS_min, E_DOS_max, dE0, E0
    INTEGER*4 :: DOS_steps

    E_DOS_MIN = -1200 * meV2au
    E_DOS_max = 1200 * meV2au
    dE0 = 1 * meV2au
    DOS_steps = INT((E_DOS_max - E_DOS_min) / dE0)


    !PRINT*, "READING INPUT"
    CALL GET_INPUT("dispersion_input.nml")

    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE(Energies(0:k1_steps, 0:k2_steps, DIM_POSITIVE_K))
    ALLOCATE(Gamma_SC(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))
    ALLOCATE(DOS(0:DOS_steps))

    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0. , 0.)
    Energies(:,:,:) = 0.
    Gamma_SC(:,:,:,:) = DCMPLX(0. , 0.)*meV2au
    Charge_dens(:) = 0.5


    !Get self consistent gamma and charge density
    !CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), "OutputData/Gamma_SC_iter.dat")
    !CALL GET_CHARGE_DENS(Charge_dens(:), "OutputData/Chargen_dens_iter.dat")


    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle


    DO i = 0, k1_steps
        DO j = 0, k2_steps
            k1 = i*dk1
            k2 = j*dk2

            kx = ( k1*SQRT(3.)/2. ) * A_TILDE
            ky = ( -k1/2. + k2 ) * A_TILDE
            Hamiltonian(:,:) = DCMPLX(0. , 0.)
            CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
            CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
            CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))
        
            CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
        
            Hamiltonian(:,:) = Hamiltonian_const(:,:) + Hamiltonian(:,:) !Should by multiplied by 0.5 if in Nambu space
        
            CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:DIM_POSITIVE_K,:DIM_POSITIVE_K), Energies(i,j,:), U_transformation(:,:), DIM_POSITIVE_K)

        END DO
    END DO

    DO n = 0, DOS_steps
        DO i = 0, k1_steps
            DO j = 0, k2_steps
                DO k = 1, DIM_POSITIVE_K
                    E0 = E_DOS_min + n*dE0
                    DOS(n) = DOS(n) + dirac_delta(Energies(i,j,k), E0)
                END DO
            END DO
        END DO
    END DO

    OPEN(unit = 9, FILE= "./OutputData/DOS.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 0, DOS_steps
        E0 = E_DOS_min + n*dE0
        WRITE(9,*) E0/meV2au, DOS(n)
    END DO
    CLOSE(9)


    CALL PRINT_ENERGIES(Energies(:,:,:), k1_steps, k2_steps, dk1, dk2, "Energies", DIM_POSITIVE_K)





    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Energies)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Charge_dens)
    DEALLOCATE(DOS)

END PROGRAM dispersion