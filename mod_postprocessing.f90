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
    !TODO: This should be specified as an input to this subroutine
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
    
    !TODO: This should be specified as an input to this subroutine
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




SUBROUTINE CALCULATE_DISPERSION(filename, brilouinZoneFraction)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER*4, INTENT(IN) :: brilouinZoneFraction
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

    !TODO: This should be specified as an input to this subroutine
    CALL GET_INPUT("dispersion_input.nml")

    yz_contribution = 0.
    zx_contribution = 0.
    xy_contribution = 0.
    lat1_contribution = 0.
    lat2_contribution = 0.
    spin_up_contribution = 0.
    spin_down_contribution = 0.


    kx_steps = INT(k1_steps/2 / brilouinZoneFraction)
    ky_steps = INT(k2_steps/2 / brilouinZoneFraction)
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

    !TODO: This should be specified as an input to this subroutine
    CALL GET_CHARGE_DENS(Charge_dens(:), "/home/jczarnecki/LAO-STO-results/LAO-STO-Hub/RUN_E_Fermi_990.0_J_SC_150.0/OutputData/Charge_dens_final.dat")

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

SUBROUTINE CALCULATE_CHERN_PARAMS(Nk1, Nk2, inputPath)
    !! Calculates Chern Params, based on https://arxiv.org/abs/cond-mat/0503172
    INTEGER*4, INTENT(IN) :: Nk1 !! Number of divisions along k1
    INTEGER*4, INTENT(IN) :: Nk2 !! Number of divisions along k2
    !INTEGER*4, INTENT(IN) :: HamDim !! DImension of the hamiltonian to be diagonalized (e.g. 4 for simple hellical, 24 for LAO-STO)
    CHARACTER(LEN=*) :: inputPath

    COMPLEX*16, ALLOCATABLE :: Psi(:,:,:,:)
    COMPLEX*16 :: U1_chern(DIM/2, DIM/2), U2_chern(DIM/2, DIM/2), U3_chern(DIM/2, DIM/2), U4_chern(DIM/2, DIM/2)
    INTEGER*4 :: i,j,a,b,m,n
    REAL*8 :: potChem
    REAL*8 :: Bfield(3)
    COMPLEX*16 :: links

    COMPLEX*16 :: f_12, det1, det2, det3, det4

    ALLOCATE(Psi(-Nk1/2:Nk1/2, 2, DIM, DIM))
    i = 0
    j = 0
    a = 0
    b = 0
    n = 0
    m = 0

    potChem = 0*meV2au !Only for testing in hellical gap
    Bfield = (/0.0*T2au, 0.0*T2au, 0.0*T2au/)

    !PRINT*, "Entered chern params"
    !PRINT*, "Ham dim ", DIM
    !PRINT*, "DIM/2", DIM/2
    !PRINT*, "Nk1/2", Nk1/2
    f_12 = DCMPLX(0., 0.)
    !Calculate Chern numbers
    DO j = -Nk2/2, Nk2/2 - 1
        !This is for memory optimization. I dont have to keep all values of Psi over Brillouin Zone.
        !Instead I need values for current row and one row above:
        ! j = 0
        ! ********************************
        ! ********************************
        ! ********************************
        ! ********************************
        ! ********************************
        ! ********************************
        ! ******************************** <- this too
        ! ******************************** <- this I need
        !In next iteration I can forget about bottom row, the second one becomes the lower one
        ! and I have to calculate one above
        ! j = 1
        ! ********************************
        ! ********************************
        ! ********************************
        ! ********************************
        ! ********************************
        ! ******************************** <- this has to be calculated, Psi (:,2,:,:)
        ! ******************************** <- this becomes Psi(:,1,:,:)
        ! ******************************** <- this I can forget
        !It could be improved to keep Nk1 + 1 values, but for now I hope it is not necessary
        IF (j .EQ. 0) THEN
            DO n = -Nk1/2, Nk1/2
                CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2, n, j, inputPath, Psi(n, 1, :, :)) !First row
                CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2, n, j+1, inputPath, Psi(n, 2, :, :)) !Second row

                !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j, Psi(n,1,:,:))
                !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j + 1, Psi(n,2,:,:))
            END DO
        ELSE
            Psi(:,1,:,:) = Psi(:,2,:,:)
            DO n = -Nk1/2, Nk1/2
                CALL LAO_STO_CHERN_ENERGIES(Nk1, Nk2,n,j+1, inputPath, Psi(n,2,:,:)) !Next row
                !CALL HELLICAL_TEST_CHERN(potChem, Bfield, Nk1, Nk2, n , j+1, Psi(n,2,:,:)) !Next row    
            END DO
        END IF

        DO i = -Nk1/2, Nk1/2 - 1
            !PRINT*, i, j
            !Calculate U matrices for chern numbers
            U1_chern = DCMPLX(0.0, 0.)
            U2_chern = DCMPLX(0.0, 0.)
            U3_chern = DCMPLX(0.0, 0.)
            U4_chern = DCMPLX(0.0, 0.)

            DO a = 1, DIM/2
                DO b = 1, DIM/2
                    U1_chern(a, b) = SUM(CONJG(Psi(i,1,:,a))*Psi(i+1,1,:,b))
                    U2_chern(a, b) = SUM(CONJG(Psi(i+1,1,:,a))*Psi(i+1,2,:,b))
                    U3_chern(a, b) = SUM(CONJG(Psi(i+1,2,:,a))*Psi(i,2,:,b))
                    U4_chern(a, b) = SUM(CONJG(Psi(i,2,:,a))*Psi(i,1,:,b))
                END DO
            END DO

            det1 = det(U1_chern(:, :), DIM/2)
            IF (det1 .ne. 0.) THEN
                det1 = det1 / ABS(det1)
            END IF

            det2 = det(U2_chern(:, :), DIM/2)
            IF (det2 .ne. 0.) THEN
                det2 = det2 / ABS(det2)
            END IF


            det3 = det(U3_chern(:, :), DIM/2)
            IF (det3 .ne. 0.) THEN
                det3 = det3 / ABS(det3)
            END IF


            det4 = det(U4_chern(:, :), DIM/2)
            IF (det4 .ne. 0.) THEN
                det4 = det4 / ABS(det4)
            END IF

            links = det1*det2*det3*det4
            f_12 = f_12 + ATAN(AIMAG(links), REAL(links))

        END DO
    END DO

    OPEN(unit = 9, FILE= TRIM(inputPath)//"OutputData/ChernNumber.dat", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,*) f_12/(2*PI)
    CLOSE(9)
    !PRINT*, "Chern number is ", f_12/(2*PI)

    DEALLOCATE(Psi)

END SUBROUTINE CALCULATE_CHERN_PARAMS

SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP(inputPath, dE, nBrillouinPoints)
    
    CHARACTER(LEN=*), INTENT(IN) :: inputPath !! This should be a path to folder where input.nml resides
    REAL*8, INTENT(IN) :: dE
    INTEGER*4, INTENT(IN) :: nBrillouinPoints
    CHARACTER(LEN=20) :: output_format
 
    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:)
    REAL*8, ALLOCATABLE :: Energies(:)

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    INTEGER*1, ALLOCATABLE :: IsFermiSurface(:,:) !! This indicates whether given (kx,ky) point is at Fermi surface
    INTEGER*2, ALLOCATABLE :: OrbitalAtFermiSurface(:,:) !! This indicates which state consitutes to the Fermi surface.
                                                         !! Because ZHEEV is used, energies are sorted from lowest to highest.
                                                         !! This implies that we do not recognize which spin, orbital etc.
                                                         !! constitutes to the Fermi surface, only "lowest", "second lowest" and so on.

    REAL*8 :: brillouinZoneVertices(6,2)

    REAL*8 :: k1, k2, kx, ky, dkx, dky
    INTEGER*4 :: i,j,n, m
    INTEGER*4 :: kx_steps, ky_steps

    LOGICAL :: fileExists

    brillouinZoneVertices(:,1) = (/ 4.*PI/(3*SQRT(3.0d0)), 2.*PI/(3*SQRT(3.0d0)), -2.*PI/(3*SQRT(3.0d0)), -4.*PI/(3*SQRT(3.0d0)), -2.*PI/(3*SQRT(3.0d0)), 2.*PI/(3*SQRT(3.0d0))/)
    brillouinZoneVertices(:,2) = (/ 0.0d0, -2.*PI/3.0d0, -2.*PI/3.0d0, 0.0d0, 2.*PI/3.0d0, 2.*PI/3.0d0/)

    dkx = KX_MAX / nBrillouinPoints
    dky = KY_MAX / nBrillouinPoints

    kx_steps = INT(nBrillouinPoints)
    ky_steps = INT(nBrillouinPoints)




    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(Energies(DIM))
    ALLOCATE(IsFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
    ALLOCATE(OrbitalAtFermiSurface(-kx_steps:kx_steps, -ky_steps:ky_steps))
    ALLOCATE(Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))


    Hamiltonian(:,:) = DCMPLX(0., 0.)
    Hamiltonian_const(:,:) = DCMPLX(0. , 0.)
    Energies(:) = 0.
    Gamma_SC(:,:,:,:) = DCMPLX(0. , 0.)*meV2au
    Charge_dens(:) = 0.
    IsFermiSurface(:,:) = 0
    OrbitalAtFermiSurface(:,:) = 0
    output_format = '(2E15.5, 2I10)'


    !Get parameters from simulation
    CALL GET_INPUT(TRIM(inputPath)//"input.nml")

    INQUIRE(FILE = TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST = fileExists)
    IF (fileExists) THEN
        CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
    ELSE
        CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
    END IF
    !Computing k-independent terms
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle

    OPEN(unit = 9, FILE= TRIM(inputPath)//"OutputData/FermiSurface.dat", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,*) '#kx[1/a] ky[1/a] is_fermi_surface N_orbital'
    !Dispersion relation in a normal state
    DO i = -kx_steps, kx_steps
        DO j = -ky_steps, ky_steps
            kx = i*dkx !* (2. * PI * 2./3.)
            ky = j*dky !* (2. * PI * 2./3.)

            IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

                Hamiltonian(:,:) = DCMPLX(0. , 0.)
                Energies(:) = 0.
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
                !CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))
            
                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
            
                Hamiltonian(:,:) = Hamiltonian_const(:,:) + Hamiltonian(:,:) !Should by multiplied by 0.5 if in Nambu space
            
                CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:DIM_POSITIVE_K, :DIM_POSITIVE_K), Energies(:DIM_POSITIVE_K), DIM_POSITIVE_K)
                !CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:DIM_POSITIVE_K, :DIM_POSITIVE_K), Energies(:DIM_POSITIVE_K), U_transformation(:DIM_POSITIVE_K, :DIM_POSITIVE_K), DIM_POSITIVE_K)

                !Check whether current wavevector is in the Fermi surface
                IF (MINVAL(ABS(Energies(:DIM_POSITIVE_K))) < dE) THEN
                    IsFermiSurface(i,j) = 1
                    OrbitalAtFermiSurface(i,j) = MINLOC(ABS(Energies(:DIM_POSITIVE_K)), 1)
                END IF
                WRITE(9, output_format) kx, ky, IsFermiSurface(i,j), OrbitalAtFermiSurface(i,j)
            END IF
        END DO
    END DO
    CLOSE(9)

    PRINT*, "Fermi surface done"

    !Calculation of superconducting gap
    INQUIRE(FILE = TRIM(inputPath)//"OutputData/Gamma_SC_final.dat", EXIST = fileExists)
    IF (fileExists) THEN
        CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
    ELSE
        CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
    END IF

    OPEN(unit = 9, FILE= TRIM(inputPath)//"OutputData/SuperconductingGap.dat", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,*) '#kx[1/a] ky[1/a] gap_SC[meV] N_orbital'
    output_format = '(3E15.5, I10)'

    DO i = -kx_steps, kx_steps
        DO j = -ky_steps, ky_steps
            IF (IsFermiSurface(i,j) == 1) THEN
                kx = i*dkx !* (2. * PI * 2./3.)
                ky = j*dky !* (2. * PI * 2./3.)

                Hamiltonian(:,:) = DCMPLX(0. , 0.)
                Energies(:) = 0.
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))
                ! DO n = 1, DIM_POSITIVE_K
                !     IF (n .le. DIM_POSITIVE_K - ORBITALS) THEN
                !         Hamiltonian(n, ORBITALS + DIM_POSITIVE_K + n) = 20 * meV2au
                !     END IF
                !     IF (n .ge. ORBITALS) THEN
                !         Hamiltonian(n, DIM_POSITIVE_K + n - ORBITALS) = 20 * meV2au
                !     END IF
                ! END DO
                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle
            
                Hamiltonian(:,:) = 0.5*(Hamiltonian_const(:,:) + Hamiltonian(:,:)) !Should by multiplied by 0.5 if in Nambu space

                CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:,:), Energies(:), DIM)
                !Write superconducting gap
                WRITE(9,output_format) kx, ky, ABS(Energies(DIM_POSITIVE_K) - Energies(DIM_POSITIVE_K + 1)) / meV2au, OrbitalAtFermiSurface(i,j)
            END IF
        END DO
    END DO
    CLOSE(9)


END SUBROUTINE CALCULATE_SUPERCONDUCTING_GAP


SUBROUTINE TRANSFORM_DELTA_MATRIX(inputPath, nBrillouinPoints)
    CHARACTER(LEN=*), INTENT(IN) :: inputPath !! This should be a path to folder where input.nml resides
    INTEGER*4, INTENT(IN) :: nBrillouinPoints

    CHARACTER(LEN=20) :: output_format
 
    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:), Hamiltonian_const(:,:), Gamma_matrix(:,:), Gamma_matrix_temp(:,:), Gamma_matrix_diagonal(:,:)
    COMPLEX*16, ALLOCATABLE :: U_transformation(:,:)
    REAL*8, ALLOCATABLE :: Energies(:)

    COMPLEX*16, ALLOCATABLE :: Gamma_SC(:,:,:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)
    REAL*8 :: brillouinZoneVertices(6,2)

    REAL*8 :: k1, k2, kx, ky, dkx, dky
    INTEGER*4 :: i,j,n, m
    INTEGER*4 :: kx_steps, ky_steps

    LOGICAL :: fileExists

    brillouinZoneVertices(:,1) = (/ 4.*PI/(3*SQRT(3.0d0)), 2.*PI/(3*SQRT(3.0d0)), -2.*PI/(3*SQRT(3.0d0)), -4.*PI/(3*SQRT(3.0d0)), -2.*PI/(3*SQRT(3.0d0)), 2.*PI/(3*SQRT(3.0d0))/)
    brillouinZoneVertices(:,2) = (/ 0.0d0, -2.*PI/3.0d0, -2.*PI/3.0d0, 0.0d0, 2.*PI/3.0d0, 2.*PI/3.0d0/)

    dkx = KX_MAX / nBrillouinPoints
    dky = KY_MAX / nBrillouinPoints

    kx_steps = INT(nBrillouinPoints)
    ky_steps = INT(nBrillouinPoints)

    ALLOCATE(Hamiltonian(DIM,DIM)) 
    ALLOCATE(Hamiltonian_const(DIM,DIM))
    ALLOCATE(U_transformation(DIM,DIM))
    ALLOCATE(Gamma_matrix(DIM,DIM))
    ALLOCATE(Gamma_matrix_diagonal(DIM,DIM))
    ALLOCATE(Gamma_matrix_temp(DIM,DIM))
    ALLOCATE(Energies(DIM))
    ALLOCATE(Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS,2, SUBLATTICES))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))

    Gamma_SC = DCMPLX(0.0d0,0.0d0)
    Charge_dens = 0.0d0

    !Get parameters from simulation
    CALL GET_INPUT(TRIM(inputPath)//"input.nml")
    !Charge densities
    ! INQUIRE(FILE = TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST = fileExists)
    ! IF (fileExists) THEN
    !     CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
    ! ELSE
    !     CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
    ! END IF
    ! !Superconducting gap parameters
    ! INQUIRE(FILE = TRIM(inputPath)//"OutputData/Gamma_SC__final.dat", EXIST = fileExists)
    ! IF (fileExists) THEN
    !     CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
    ! ELSE
    !     CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
    ! END IF

    Gamma_SC(:,:,1,:) = 999. * meV2au
    Gamma_SC(:,:,2,:) = -999. * meV2au

    !Initialize constant hamiltonian
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
            kx = i*dkx !* (2. * PI * 2./3.)
            ky = j*dky !* (2. * PI * 2./3.)
            
            IF (is_inside_polygon(brillouinZoneVertices, 6, kx, ky)) THEN

                Hamiltonian(:,:) = DCMPLX(0. , 0.)
                Energies(:) = 0.
                CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), kx, ky)
                CALL COMPUTE_TI1_TI2(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_PI(Hamiltonian(:,:), kx, ky) !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_H_SIGMA(Hamiltonian(:,:), kx, ky)  !There may be a problem since Ti1,Ti2 coupling is assumed to be equal Ti2,Ti1
                CALL COMPUTE_HUBBARD(Hamiltonian(:,:), Charge_dens(:))
                CALL COMPUTE_SC(Hamiltonian(:,:), kx, ky, Gamma_SC(:,:,:,:))
                CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle

                Hamiltonian(:,:) = 0.5*(Hamiltonian_const(:,:) + Hamiltonian(:,:)) !Should by multiplied by 0.5 if in Nambu space
                ! DO n = 1, DIM
                !     WRITE(*,'(24E12.4)') (ABS(Hamiltonian(n,m))**2/meV2au, m = 1, DIM)
                ! END DO
                ! WRITE(*,*)
                ! WRITE(*,*)

                !CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:, :), Energies(:), DIM)
                CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), DIM)

                !Compute current superconducing coupling matrix separately
                Gamma_matrix(:,:) = DCMPLX(0.0d0, 0.0d0)
                Gamma_matrix_diagonal(:,:) = DCMPLX(0.0d0, 0.0d0)
                Gamma_matrix_temp(:,:) = DCMPLX(0.0d0, 0.0d0)
                CALL COMPUTE_SC(Gamma_matrix(:,:), kx, ky, Gamma_SC(:,:,:,:))
                CALL COMPUTE_CONJUGATE_ELEMENTS(Gamma_matrix(:,:), DIM)
                Gamma_matrix = 0.5 * Gamma_matrix
                DO n = 1, DIM
                    DO m = 1, DIM
                        Gamma_matrix_temp(n,m) = SUM(Gamma_matrix(n,:)*U_transformation(:,m))
                    END DO
                END DO
                
                DO n = 1, DIM
                    DO m = 1, DIM
                        Gamma_matrix_diagonal(n,m) = SUM(CONJG(U_transformation(:,n))*Gamma_matrix_temp(:,m))
                    END DO
                END DO

                !Gamma_matrix_diagonal = MATMUL(MATMUL(Gamma_matrix, U_transformation), TRANSPOSE(CONJG(U_transformation)))

                WRITE(*,'(A, 24F10.4)') 'Energies: ', (Energies(m)/meV2au, m = 1, DIM)
                DO n = 1, DIM
                    WRITE(*,'(24F10.4)') (REAL(Gamma_matrix_diagonal(n,m))/meV2au, m = 1, DIM)
                END DO
                WRITE(*,*)
                WRITE(*,*)

            END IF
        END DO
    END DO

END SUBROUTINE TRANSFORM_DELTA_MATRIX

SUBROUTINE HELLICAL_TEST_CHERN(potChem, B, Nk1, Nk2, i, j, U_transformation)
    INTEGER*4, PARAMETER :: HamDim = 4
    REAL*8, PARAMETER :: g = 5.0d0
    REAL*8, PARAMETER :: tHop = 200*meV2au
    REAL*8, PARAMETER :: alphaSOC = 100*meV2au
    REAL*8, PARAMETER :: gammaSC = DCMPLX(0.5*meV2au, 0.0d0)
    REAL*8, PARAMETER :: muB = 0.5d0
    
    
    REAL*8, INTENT(IN) :: potChem
    REAL*8, INTENT(IN) :: B(3) ![Bx, By, Bz]
    INTEGER*4, INTENT(IN) :: Nk1, Nk2, i ,j

    COMPLEX*16, INTENT(INOUT) :: U_transformation(HamDim, HamDim)
    REAL*8 :: Energies(HamDim)
    REAL*8 :: dkx, dky
    REAL*8 :: kx, ky
    INTEGER*4 :: m,n

    COMPLEX*16 :: Hamiltonian(HamDim, HamDim)

    dkx = 2.0d0*PI/Nk1
    dky = 2.0d0*PI/Nk2

    Energies(:) = 0.0d0
    U_transformation(:,:) = DCMPLX(0.0d0, 0.0d0)
    Hamiltonian(:,:) = 0.0d0
    !Calculate all points in Brillouin zone
    kx = i*dkx
    ky = j*dky

    !Diagonal terms
    !Electrons H(k)
    Hamiltonian(1, 1) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem + 0.5*muB*g*B(3)
    Hamiltonian(2, 2) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem - 0.5*muB*g*B(3)
    
    !Holes -H*(-k)
    Hamiltonian(3, 3) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem + 0.5*muB*g*B(3))
    Hamiltonian(4, 4) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem - 0.5*muB*g*B(3))

    !Spin-orbit coupling
    Hamiltonian(1,2) = 0.5*mub*g*(B(1) - imag*B(2)) + alphaSOC*(DSIN(kx) + imag*DSIN(ky))
    Hamiltonian(3,4) = -(0.5*mub*g*(B(1) + imag*B(2)) + alphaSOC*(DSIN(-kx) - imag*DSIN(-ky)))

    !Superconductivity
    Hamiltonian(1, 4) = gammaSC
    Hamiltonian(2,3) = -gammaSC

    Hamiltonian(:,:) = 0.5 * Hamiltonian(:,:)
    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, HamDim)

    CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), HamDim)

    CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies(:), U_transformation(:,:), HamDim)

END SUBROUTINE HELLICAL_TEST_CHERN

SUBROUTINE LAO_STO_CHERN_ENERGIES(Nk1, Nk2, i, j, inputPath, U_transformation)
    !! This subroutine calculates energies and wavefunctions of LAO-STO in [111] direction.
    !! Returns sorted wavefunctions in (i,j) point of the Brillouin zone.
    INTEGER*4, INTENT(IN) :: Nk1 !! Number of divisions of Brillouin zone in direction k1.
    INTEGER*4, INTENT(IN) :: Nk2 !! Number of divisions of Brillouin zone in direction k2.
    INTEGER*4, INTENT(IN) :: i !! Curent point k1_i in the Brillouin zone
    INTEGER*4, INTENT(IN) :: j !! Curent point k2_j in the Brillouin zone
    CHARACTER(LEN=*), INTENT(IN) :: inputPath !! Directory of the run, where input.nml should be placed.
                                            !! It contains material information and physical parameter of calculation:
                                            !! Fermi energy, temperature etc.
    COMPLEX*16, INTENT(OUT) :: U_transformation(DIM, DIM) !! Matrix containing eigenvectors stored in consecutive columns.
                                                          !! On output sorted based on energies from lowest to highest.
    
    REAL*8 :: Energies(DIM) !! Eigenvalues of the hamiltonian

    COMPLEX*16 :: Hamiltonian(DIM,DIM), Hamiltonian_const(DIM,DIM)
    COMPLEX*16 :: Gamma_SC(ORBITALS, N_ALL_NEIGHBOURS, 2, SUBLATTICES)
    REAL*8 :: Charge_dens(DIM_POSITIVE_K)
    REAL*8 :: k1, k2, kx, ky
    REAL*8 :: dk1_Chern, dk2_Chern
    INTEGER*4 :: n
    LOGICAL :: fileExists

    !PRINT*, "Allocation ended"
    dk1_Chern = K1_MAX/Nk1
    dk2_Chern = K2_MAX/Nk2

    U_transformation = DCMPLX(0. , 0.)
    Energies = 0.
    Gamma_SC = 0.
    Charge_dens = 0.

    !PRINT*, "Entered chern energies"
    !Get parameters from simulation
    CALL GET_INPUT(TRIM(inputPath)//"input.nml")
    !CHeck if units are correct
    !Get Gammas from run
    INQUIRE(FILE = TRIM(inputPath)//"OutputData/Gamma_SC__final.dat", EXIST = fileExists)
    IF (fileExists) THEN
        CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_final.dat")
    ELSE
        CALL GET_GAMMA_SC(Gamma_SC(:,:,:,:), TRIM(inputPath)//"OutputData/Gamma_SC_iter.dat")
    END IF
    
    !Get cherge densities from file
    INQUIRE(FILE = TRIM(inputPath)//"OutputData/Charge_dens_final.dat", EXIST = fileExists)
    IF (fileExists) THEN
        CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_final.dat")
    ELSE
        CALL GET_CHARGE_DENS(Charge_dens(:), TRIM(inputPath)//"OutputData/Charge_dens_iter.dat")
    END IF

    !Gamma_SC(:,:,1,:) = 100.0d0 * meV2au
    !Gamma_SC(:,:,2,:) = -100.0d0 * meV2au


    !Computing k-independent terms
    Hamiltonian_const = DCMPLX(0., 0.)
    CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian_const(:,:))
    CALL COMPUTE_ELECTRIC_FIELD(Hamiltonian_const(:,:))
    !CALL COMPUTE_ZEEMAN((/0.0d0, 0.0d0, 0.1*T2au/), Hamiltonian_const(:,:))
    DO n = 1, DIM_POSITIVE_K
        Hamiltonian_const(n,n) = Hamiltonian_const(n,n) - E_Fermi
        Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) = Hamiltonian_const(DIM_POSITIVE_K + n, DIM_POSITIVE_K + n) + E_Fermi
    END DO

    CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian_const(:,:), DIM) !This is not needed, since ZHEEV takes only upper triangle

    !Calculate eigenvalues and eigenvectors to later compute chern numbers
    k1 = i*dk1_Chern
    k2 = j*dk2_Chern

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

    CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), DIM)

    CALL SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, U_transformation, DIM)

    ! PRINT*, "ENERGIES FOR ", i, j, k1, k2
    ! DO n = 1, DIM
    !     PRINT*, Energies(n)
    ! END DO

END SUBROUTINE LAO_STO_CHERN_ENERGIES


!########################### HELPER FUNCTIONS ###########################
SUBROUTINE SORT_ENERGIES_AND_WAVEFUNCTIONS(Energies, Psi, HamDim)
    INTEGER*4, INTENT(IN) :: HamDim
    COMPLEX*16, INTENT(INOUT) :: Psi(HamDim, HamDim)
    REAL*8, INTENT(INOUT) :: Energies(HamDim)

    INTEGER*4 :: i,j
    REAL*8 :: tmpEnergy
    COMPLEX*16 :: tmpPsi(HamDim)

    DO i = 1, HamDim
        DO j = 1, HamDim - 1
        IF (Energies(j) .GT. Energies(j + 1)) THEN
            !Swap energies
            tmpEnergy = Energies(j)
            Energies(j) = Energies(j + 1)
            Energies(j + 1) = tmpEnergy

            !Swap wavefunctions
            tmpPsi(:) = Psi(:, j)
            Psi(:, j) = Psi(:, j + 1)
            Psi(:, j + 1) = tmpPsi
        END IF
        END DO
    END DO

    DO i = 1, HamDim
        Psi(:,i) = Psi(:,i)/SUM(ABS(Psi(:,i))**2)
    END DO

    ! DO i = 1, HamDim
    !     !PRINT*, "Psi(1,i) ", Psi(1,i)
    !     Psi(:,i) = Psi(:,i) / (Psi(1,i) / ABS(Psi(i,i)))
    ! END DO

END SUBROUTINE SORT_ENERGIES_AND_WAVEFUNCTIONS

COMPLEX*16 FUNCTION det(matrix, n)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: n
    COMPLEX*16, INTENT(IN) :: matrix(n,n)
    INTEGER*4 :: IPIV(n)
    INTEGER*4 :: info, i
    
    IPIV(:) = 0.0d0

    CALL ZGETRF(n,n,matrix,n,IPIV, info)
    CALL ZLAPMT(.TRUE., n, n, matrix, n, IPIV)
    det = 1.0d0
    DO i = 1, n
        det = det * matrix(i,i)
    END DO
    RETURN

END FUNCTION det

LOGICAL FUNCTION is_inside_polygon(verticesArray, nVertices, pointX, pointY)
    !! This function returns true if the point (pointX, pointY) is inside the polygon
    !! Where vertices are defined in array verticesArray, containing X and Y coordinates.
    !! The assumption is that the polygon is two-dimensional
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nVertices !! Number of vertices of polygon
    REAL*8, INTENT(IN) :: verticesArray(nVertices,2) !! X and Y coordinates of vertices
    REAL*8, INTENT(IN) :: pointX, pointY !! X and Y coordinates of point to be tested
    INTEGER*4 :: i, j
    REAL*8 :: xIntersection

    is_inside_polygon = .FALSE.

    DO i = 1, nVertices

        IF (pointX == verticesArray(i,1) .AND. pointY == verticesArray(i,2)) THEN
            is_inside_polygon = .TRUE.
            RETURN
        END IF

        j = MOD(i, nVertices) + 1 !! Index of next vertex
        !! Check if horizontal line of y = pointY can cross the line joining two vertices
        IF (pointY > MIN(verticesArray(i,2), verticesArray(j,2)) .AND. &
        & pointY <= MAX(verticesArray(i,2), verticesArray(j,2))) THEN
            IF (pointX <= MAX(verticesArray(i,1), verticesArray(j,1)) ) THEN
                IF (verticesArray(i,2) /= verticesArray(j,2)) THEN !To avoid division by zero
                    !Calculates intersection point of line connecting two vertices and horizontal line y = pointY
                    xIntersection = (pointY - verticesArray(i,2)) * &
                    & (verticesArray(j,1) - verticesArray(i,1))/(verticesArray(j,2) &
                    & - verticesArray(i,2)) + verticesArray(i,1)

                    IF (pointX <= xIntersection .OR. verticesArray(i,1) == verticesArray(j,1)) THEN
                        is_inside_polygon = .NOT. is_inside_polygon
                    END IF
                END IF
            END IF
        END IF
    END DO

    RETURN
END FUNCTION is_inside_polygon

    !MAKE THIS A SUBROUTINE
    !Plot this dispersion
    ! OPEN(unit = 9, FILE= "./OutputData/EnergiesHellical.dat", FORM = "FORMATTED", ACTION = "WRITE")
    ! DO n = -Nk1/2, Nk1/2
    !     kx = n*dkx
    !     ky = kx

    !     ! Hamiltonian(1, 1) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) + 0.5*muB*g*B(3)
    !     ! Hamiltonian(2, 2) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - 0.5*muB*g*B(3)
        
    !     ! !Spin-orbit coupling
    !     ! Hamiltonian(1,2) = 0.5d0*mub*g*(B(1) - imag*B(2)) + alphaSOC*(DSIN(kx) + imag*DSIN(ky))
    !     ! CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian(:2, :2), 2)

    !     ! CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:2, :2), Energies(:2), U_transformation(:2,:2), 2)
        
    !     !Diagonal terms
    !     !Electrons H(k)
    !     Hamiltonian(1, 1) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem + 0.5*muB*g*B(3)
    !     Hamiltonian(2, 2) = 2.0*tHop*(1. - DCOS(kx)) + 2.0*tHop*(1 - DCOS(ky)) - potChem - 0.5*muB*g*B(3)
        
    !     !Holes -H*(-k)
    !     Hamiltonian(3, 3) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem + 0.5*muB*g*B(3))
    !     Hamiltonian(4, 4) = -(2.0*tHop*(1. - DCOS(-kx)) + 2.0*tHop*(1 - DCOS(-ky)) - potChem - 0.5*muB*g*B(3))
    
    !     !Spin-orbit coupling
    !     Hamiltonian(1,2) = 0.5*mub*g*(B(1) - imag*B(2)) + alphaSOC*(DSIN(kx) + imag*DSIN(ky))
    !     Hamiltonian(3,4) = -(0.5*mub*g*(B(1) + imag*B(2)) + alphaSOC*(DSIN(-kx) - imag*DSIN(-ky)))
    
    !     !Superconductivity
    !     Hamiltonian(1, 4) = gammaSC
    !     Hamiltonian(2,3) = -gammaSC
    
    !     Hamiltonian(:,:) = 0.5 * Hamiltonian(:,:)
    !     CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, HamDim)
    
    !     CALL DIAGONALIZE_GENERALIZED(Hamiltonian(:,:), Energies(:), U_transformation(:,:), HamDim)
   
    !     DO m = 1, HamDim
    !         WRITE(9,*) kx, Energies(m)/meV2au
    !     END DO

    ! END DO
    ! CLOSE(9)

END MODULE mod_postprocessing