PROGRAM chern
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_compute_hamiltonians
    IMPLICIT NONE
    
    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:,:,:)
    COMPLEX*16, ALLOCATABLE :: U1_chern(:,:), U2_chern(:,:), U3_chern(:,:), U4_chern(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)
    REAL*8 :: k1, k2, kx, ky
    INTEGER*4 :: i,j,k,n,a,b, lat, orb, orb_prime,spin, a_prime, b_prime
    CHARACTER(LEN=100) run_dir
    REAL*8 :: step_multiplier

    COMPLEX*16 :: f_12, det1, det2, det3, det4
    INTEGER*4, ALLOCATABLE :: Below_fermi(:)
    INTEGER*4 :: occupied

    COMPLEX*16, EXTERNAL :: det

    run_dir = "~/LAO-STO-results/RUNS_low_U/RUN_E_Fermi_-800.0_U_HUB_0.0_V_HUB_0.0/"
    CALL GET_INPUT(TRIM(run_dir)//"input.nml")

    !Modify k_steps and dk to refine grid for Chern number calculations
    step_multiplier = 20
    k1_steps = k1_steps*step_multiplier
    k2_steps = k2_steps*step_multiplier
    dk1 = dk1/step_multiplier
    dk2 = dk2/step_multiplier

    ALLOCATE(Hamiltonian(DIM,DIM))
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(0:k1_steps,0:k2_steps, DIM, DIM))
    ALLOCATE(U1_chern(DIM,DIM))
    ALLOCATE(U2_chern(DIM,DIM))
    ALLOCATE(U3_chern(DIM,DIM))
    ALLOCATE(U4_chern(DIM,DIM))
    ALLOCATE(Below_fermi(DIM))
    ALLOCATE(Energies(0:k1_steps,0:k2_steps, DIM))
    ALLOCATE(Gamma_SC(ORBITALS,N_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))


    U_transformation = DCMPLX(0. , 0.)
    Energies = 0.
    Gamma_SC = 0.
    Charge_dens = 0.

    CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(run_dir)//"OutputData/Gamma_SC_final.dat")
    CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(run_dir)//"OutputData/Chargen_dens_final.dat")

    Gamma_SC = Gamma_SC*100

    !Computing k-independent terms
    Hamiltonian_const = DCMPLX(0., 0.)
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
  
    !Calculate eigenvalues and eigenvectors to later compute chern numbers
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
        
            Hamiltonian(:,:) = 0.5*(Hamiltonian_const(:,:) + Hamiltonian(:,:)) !Should by multiplied by 0.5 if in Nambu space
        
            CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(i,j,:), U_transformation(i,j,:,:), DIM)

            !Normalize wavefunctions
            DO n = 1, DIM
                U_transformation(i,j,:,n) = U_transformation(i,j,:,n)/SUM(ABS(U_transformation(i,j,:,n))**2)
            END DO
        END DO
    END DO


    f_12 = DCMPLX(0., 0.)
    !Calculate Chern numbers
    DO i = 0, k1_steps - 1
        DO j = 0, k2_steps - 1
            
            !Calculate U matrices for chern numbers
            U1_chern = DCMPLX(0.0, 0.)
            U2_chern = DCMPLX(0.0, 0.)
            U3_chern = DCMPLX(0.0, 0.)
            U4_chern = DCMPLX(0.0, 0.)

            !Check what energies are below Fermi energy
            !This is an approximation, because we consider that if E(kx,ky) is below Fermi energy
            !then also E(kx + dk, ky + dk) etc. would be below Fermi energy
            Below_fermi(:) = 0
            DO n = 1, DIM
                IF(Energies(i,j,n) < 0) Below_fermi(n) = 1
            END DO

            a_prime = 0
            DO a = 1, DIM
                IF (Below_fermi(a) ) THEN
                    a_prime = a_prime + 1
                    b_prime = 0
                    DO b = 1, DIM
                        IF (Below_fermi(b)) THEN
                            b_prime = b_prime + 1
                            U1_chern(a_prime,b_prime) = SUM(CONJG(U_transformation(i,j,:,a))*U_transformation(i+1,j,:,b))
                            U2_chern(a_prime,b_prime) = SUM(CONJG(U_transformation(i+1,j,:,a))*U_transformation(i+1,j+1,:,b))
                            U3_chern(a_prime,b_prime) = SUM(CONJG(U_transformation(i,j+1,:,a))*U_transformation(i+1,j+1,:,b))
                            U4_chern(a_prime,b_prime) = SUM(CONJG(U_transformation(i,j,:,a))*U_transformation(i,j+1,:,b))        
                        END IF
                    END DO
                END IF
                !PRINT*, a_prime, b_prime, Energies(i,j,a)
            END DO

            occupied = SUM(Below_fermi(:))

            det1 = det(U1_chern(:occupied, :occupied), occupied)
            det1 = det1 / ABS(det1)

            det2 = det(U2_chern(:occupied, :occupied), occupied)
            det2 = det2 / ABS(det2)

            det3 = det(U3_chern(:occupied, :occupied), occupied)
            det3 = det3 / ABS(det3)

            det4 = det(U4_chern(:occupied, :occupied), occupied)
            det4 = det4 / ABS(det4)

            f_12 = f_12 + LOG(det1*det2/det3/det4)
        END DO
    END DO

    PRINT*, "Chern number is ", f_12*dk1*dk2/(K1_MAX*K2_MAX*imag)





    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(U_transformation)
    DEALLOCATE(U1_chern)
    DEALLOCATE(U2_chern)
    DEALLOCATE(U3_chern)
    DEALLOCATE(U4_chern)
    DEALLOCATE(Below_fermi)
    DEALLOCATE(Energies)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Charge_dens)

END PROGRAM chern

COMPLEX*16 FUNCTION det(matrix, n)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: n
    COMPLEX*16, INTENT(IN) :: matrix(n,n)
    INTEGER*4 :: IPIV(n)
    INTEGER*4 :: info, i
    
    IPIV(:) = 0. 

    CALL ZGETRF(n,n,matrix,n,IPIV, info)
    CALL ZLAPMT(.TRUE., n, n, matrix, n, IPIV)
    det = 1.
    DO i = 1, n
        det = det * matrix(i,i)
    END DO
    RETURN

END FUNCTION det