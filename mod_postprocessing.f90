MODULE mod_postprocessing
USE mod_hamiltonians
USE mod_parameters
USE mod_utilities
USE mod_writers
USE mod_reader
USE mod_compute_hamiltonians
IMPLICIT NONE
CONTAINS


SUBROUTINE CALCULATE_DOS(E_DOS_min, E_DOS_max, dE0, zeta_DOS, DOS_filename)

    !DOS calculation
    REAL*8, INTENT(IN) :: E_DOS_min, E_DOS_max, dE0
    REAL*8, INTENT(IN) :: zeta_DOS
    CHARACTER(LEN=*), INTENT(IN) :: DOS_filename
    REAL*8 :: E0
    INTEGER*4 :: DOS_steps

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    REAL*8, ALLOCATABLE :: DOS(:)
    CHARACTER(LEN=20) :: output_format


    REAL*8 :: k1, k2, kx, ky
    INTEGER*4 :: i,j,k,n, lat, orb, orb_prime,spin


    ! E_DOS_MIN = E_DOS_min * meV2au
    ! E_DOS_max = E_DOS_max * meV2au
    ! dE0 = dE0 * meV2au
    DOS_steps = INT((E_DOS_max - E_DOS_min) / dE0)
    PRINT*, E_DOS_min, E_DOS_min, dE0
    PRINT*, zeta_DOS
    PRINT*, DOS_filename

    !PRINT*, "READING INPUT"
    CALL GET_INPUT("dispersion_input.nml")

    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE(Energies(0:k1_steps, 0:k2_steps, DIM_POSITIVE_K))
    ALLOCATE(Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))
    ALLOCATE(DOS(0:DOS_steps))

    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0. , 0.)
    Energies(:,:,:) = 0.
    Gamma_SC(:,:,:,:) = DCMPLX(0. , 0.)*meV2au
    Charge_dens(:) = 0.


    !Get self consistent gamma and charge density
    !CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), "OutputData/Gamma_SC_iter.dat")
    CALL GET_CHARGE_DENS(Charge_dens(:), "/home/jczarnecki/LAO-STO-results/RUNS_low_U/RUN_E_Fermi_-905.0_U_HUB_166.66666666666666_V_HUB_166.66666666666666/OutputData/Chargen_dens_final.dat")


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

            kx = 2.*PI/(SQRT(3.0d0)) * k1
            ky = -2.*PI/3. * k1 + 4.*PI/3. * k2
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
                    DOS(n) = DOS(n) + dirac_delta(Energies(i,j,k), E0, zeta_DOS)
                END DO
            END DO
        END DO
    END DO


    output_format = '(2E15.5)'
    OPEN(unit = 9, FILE= TRIM(DOS_filename), FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,'(A)') "#E[meV] DOS[a.u]"
    DO n = 0, DOS_steps
        E0 = E_DOS_min + n*dE0
        WRITE(9,output_format) E0/meV2au, DOS(n)
    END DO
    CLOSE(9)

    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Energies)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Charge_dens)
    DEALLOCATE(DOS)

END SUBROUTINE CALCULATE_DOS




SUBROUTINE CALCULATE_DISPERSION(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)
    REAL*8, ALLOCATABLE :: Probability(:,:,:,:) ! |Psi^2|

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    REAL*8, ALLOCATABLE :: DOS(:)

    REAL*8 :: k1, k2, kx, ky
    INTEGER*4 :: i,j,k,n, lat, orb, orb_prime,spin, l, m
    INTEGER*4 :: kx_steps, ky_steps
    REAL*8 :: yz_contribution, zx_contribution, xy_contribution
    REAL*8 :: lat1_contribution, lat2_contribution
    REAL*8 :: spin_up_contribution, spin_down_contribution

    !PRINT*, "READING INPUT"
    CALL GET_INPUT("dispersion_input.nml")

    yz_contribution = 0.
    zx_contribution = 0.
    xy_contribution = 0.
    lat1_contribution = 0.
    lat2_contribution = 0.
    spin_up_contribution = 0.
    spin_down_contribution = 0.


    kx_steps = INT(k1_steps/2)
    ky_steps = INT(k2_steps/2)
    output_format = '(I5, 10E15.5)'

    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE(Probability(-kx_steps:kx_steps, -ky_steps:ky_steps, DIM_POSITIVE_K, DIM_POSITIVE_K))
    ALLOCATE(Energies(-kx_steps:kx_steps, -ky_steps:ky_steps, DIM_POSITIVE_K))
    ALLOCATE(Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))

    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0. , 0.)
    Probability(:,:,:,:) = 0.
    Energies(:,:,:) = 0.
    Gamma_SC(:,:,:,:) = DCMPLX(0. , 0.)*meV2au
    Charge_dens(:) = 0.
    CALL GET_CHARGE_DENS(Charge_dens(:), "/home/jczarnecki/LAO-STO-results/RUNS_low_U/RUN_E_Fermi_-905.0_U_HUB_166.66666666666666_V_HUB_166.66666666666666/OutputData/Chargen_dens_final.dat")

    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle

    DO i = -kx_steps, kx_steps
        DO j = -ky_steps, ky_steps
            kx = i*dk1 * (2. * PI * 2./3.)
            ky = j*dk2 * (2. * PI * 2./3.)

            Hamiltonian(:,:) = DCMPLX(0. , 0.)
            CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
            CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
            CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
            CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))
        
            CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
        
            Hamiltonian(:,:) = Hamiltonian_const(:,:) + Hamiltonian(:,:) !Should by multiplied by 0.5 if in Nambu space
        
            !CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:DIM_POSITIVE_K,:DIM_POSITIVE_K), Energies(i,j,:), U_transformation(:,:), DIM_POSITIVE_K)
            !Probability(i,j,:,:) = ABS(U_transformation)**2

            CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:DIM_POSITIVE_K, :DIM_POSITIVE_K), Energies(i,j,:), DIM_POSITIVE_K)
            Probability(i,j,:,:) = ABS(Hamiltonian(:DIM_POSITIVE_K,:DIM_POSITIVE_K))**2
            
        END DO
    END DO

    OPEN(unit = 9, FILE= "./OutputData/"//filename, FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,'(A)') "#N kx[1/a] ky[1/a] Energy[meV] P(yz) P(zx) P(xy) P(lat1) P(lat2) P(s_up) P(s_down)"
    DO l = 1, DIM_POSITIVE_K
        DO i = -kx_steps, kx_steps
            DO j = -ky_steps, ky_steps
                kx = i*dk1 * (2. * PI * 2./3.)
                ky = j*dk2 * (2. * PI * 2./3.)


                !Calculate specific contributions
                !Distinguishing orbital contributions
                yz_contribution = 0.
                zx_contribution = 0.
                xy_contribution = 0.
                DO n = 1, DIM_POSITIVE_K, 3
                    yz_contribution = yz_contribution + Probability(i,j,n,l)
                    zx_contribution = zx_contribution + Probability(i,j,n+1,l)
                    xy_contribution = xy_contribution + Probability(i,j,n+2,l)
                END DO

                !Distinguishing lattice contributions
                lat1_contribution = 0.
                lat2_contribution = 0.                
                DO n = 1, 3
                    DO m = 0, 1
                        lat1_contribution = lat1_contribution + Probability(i,j,m*TBA_DIM + n,l)
                        lat2_contribution = lat2_contribution + Probability(i,j,m*TBA_DIM + n + 3,l)
                    END DO
                END DO

                !Distinguishing spin contributions
                spin_up_contribution = 0.
                spin_down_contribution = 0.
                DO n = 1, TBA_DIM
                    spin_up_contribution = spin_up_contribution + Probability(i,j,n,l)
                    spin_down_contribution = spin_down_contribution + Probability(i,j,TBA_DIM + n,l)
                END DO

                WRITE(9, output_format) l, kx, ky, Energies(i, j, l)/meV2au, &
                & yz_contribution, zx_contribution, xy_contribution, &
                & lat1_contribution, lat2_contribution, &
                & spin_up_contribution, spin_down_contribution
            END DO
        END DO
        WRITE(9,*)
        WRITE(9,*)
    END DO 
    CLOSE(9)   





    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Hamiltonian_const)
    DEALLOCATE(U_transformation)
    DEALLOCATE(Probability)
    DEALLOCATE(Energies)
    DEALLOCATE(Gamma_SC)
    DEALLOCATE(Charge_dens)

END SUBROUTINE CALCULATE_DISPERSION



SUBROUTINE CALCULATE_CHERN()
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

    !COMPLEX*16, EXTERNAL :: determinant

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
    ALLOCATE(Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES))
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

            kx = 2.*PI/(SQRT(3.0d0)) * k1
            ky = -2.*PI/3. * k1 + 4.*PI/3. * k2
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

    PRINT*, "Chern number is ", f_12*dk1*dk2/imag

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

END SUBROUTINE CALCULATE_CHERN


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



END MODULE mod_postprocessing