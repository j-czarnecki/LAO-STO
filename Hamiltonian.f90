PROGRAM Article
    IMPLICIT NONE 
    SAVE 

    COMPLEX*16, EXTERNAL :: epsilon_yz, epsilon_zx, epsilon_xy
    REAL*8, EXTERNAL :: meV_to_au, nm_to_au, au_to_meV, au_to_nm

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:)
    REAL*8, ALLOCATABLE :: Energies(:,:,:)
    
    
    INTEGER*4, PARAMETER :: DIM = 12 !Hamiltonian is 12x12
    INTEGER*4, PARAMETER :: TBA_DIM = 6
    COMPLEX*16, PARAMETER :: imag = DCMPLX(0. , 1.)


    REAL*8 :: t_D = 1e3 * 0.5 !meV
    REAL*8 :: t_I = 1e3 * 0.04 !meV
    REAL*8:: lambda_SOC = 1e3 * 0.01 !meV
    REAL*8 :: A_TILDE =  SQRT(2./3.)*0.3905 !nm
    REAL*8 :: DELTA_TRI =  -0.005 * 1e3 !meV
    REAL*8 :: kx, ky 

    REAL*8 :: dk, k_max !Suppose kx_max == ky_max and dk_x == dk_y
    INTEGER*4 :: k_steps

    INTEGER*4 :: i,j,n

    PRINT*, A_TILDE
    !Conversion to atomic units
    t_D = meV_to_au(t_D)
    t_I = mev_to_au(t_I)
    lambda_SOC = mev_to_au(lambda_SOC)
    A_TILDE = nm_to_au( A_TILDE)
    DELTA_TRI = mev_to_au(DELTA_TRI)

    PRINT*, A_TILDE

    k_max = 10. * 1./A_TILDE !Multiplier is arbitrary. Part of Brillouin zone to integrate over
    k_steps = 2000
    dk = k_max / k_steps



    ALLOCATE(Hamiltonian(DIM,DIM))
    ALLOCATE(Energies(-k_steps:k_steps, -k_steps:k_steps, DIM))

    OPEN(unit = 9, FILE= "./E_k.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO i = -k_steps, k_steps
        DO j = 0, 0
            kx = i*dk*A_TILDE
            ky = j*dk*A_TILDE
            Energies(i,j,:) = 0.
            Hamiltonian(:,:) = 0.
            CALL COMPUTE_TBA_TERM(Hamiltonian(:,:), TBA_DIM, kx, ky, t_D, t_I)
            CALL COMPUTE_TRIGONAL_TERMS(Hamiltonian(:,:), TBA_DIM, DELTA_TRI)
            CALL COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian(:,:), DIM, lambda_SOC)

            CALL COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, DIM)
            CALL DIAGONALIZE_HERMITIAN(Hamiltonian(:,:), Energies(i,j,:), DIM)
            
            DO n = 1, DIM 
                WRITE(9,*) kx/A_TILDE, ky/A_TILDE, au_to_meV(Energies(i,j,n))*1e-3
            END DO 
            WRITE(9,*)
            WRITE(9,*)
        END DO
    END DO
    CLOSE(9)


    OPEN(unit = 9, FILE= "./E_k_sorted.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 1, DIM
        DO i = -k_steps, k_steps
            WRITE(9,*) i*dk, 0, au_to_meV(Energies(i, 0, n))*1e-3
        END DO
        WRITE(9,*)
        WRITE(9,*)
    END DO 
    CLOSE(9)


    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Energies)
    
END PROGRAM Article

SUBROUTINE COMPUTE_TBA_TERM(Hamiltonian, TBA_DIM, kx, ky, t_D, t_I)
    IMPLICIT NONE 
    INTEGER*4 :: TBA_DIM
    REAL*8 :: t_D, t_I, kx, ky
    COMPLEX*16 :: Hamiltonian(2*TBA_DIM,2*TBA_DIM) !Twice as big because of spin
    COMPLEX*16, EXTERNAL :: epsilon_yz, epsilon_zx, epsilon_xy

    !Only specifying upper triangle of matrix, since Hamiltonian is hermitian
    !spin up
    Hamiltonian(1,4) = epsilon_yz(t_D, t_I, kx, ky) 
    Hamiltonian(2,5) = epsilon_zx(t_D, t_I, kx, ky) 
    Hamiltonian(3,6) = epsilon_xy(t_D, t_I, kx, ky) 
    !spin down
    Hamiltonian(TBA_DIM + 1,TBA_DIM + 4) = epsilon_yz(t_D, t_I, kx, ky)
    Hamiltonian(TBA_DIM + 2,TBA_DIM + 5) = epsilon_zx(t_D, t_I, kx, ky)
    Hamiltonian(TBA_DIM + 3,TBA_DIM + 6) = epsilon_xy(t_D, t_I, kx, ky)
END SUBROUTINE COMPUTE_TBA_TERM

SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian, DIM, lambda_SOC)
    IMPLICIT NONE 
    INTEGER*4, INTENT(IN):: DIM
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)
    REAL*8, INTENT(IN) :: lambda_SOC
    COMPLEX*16, PARAMETER :: imag = DCMPLX(0. , 1.)

    Hamiltonian(1,2) = Hamiltonian(1,2) + lambda_SOC/2.
    Hamiltonian(1,9) = Hamiltonian(1,9) - lambda_SOC/2.
    Hamiltonian(2,9) = Hamiltonian(2,9) + imag*lambda_SOC/2.
    Hamiltonian(3,7) = Hamiltonian(3,7) + lambda_SOC/2.
    Hamiltonian(3,8) = Hamiltonian(3,8) - imag*lambda_SOC/2.
    Hamiltonian(4,5) =  Hamiltonian(4,5) + imag*lambda_SOC/2.
    Hamiltonian(4,12) = Hamiltonian(4,12) - lambda_SOC/2.
    Hamiltonian(5,12) =  Hamiltonian(5,12) + imag*lambda_SOC/2.
    Hamiltonian(6,10) =  Hamiltonian(6,10) + lambda_SOC/2.
    Hamiltonian(6,11) = Hamiltonian(6,11) - imag*lambda_SOC/2.
    Hamiltonian(7,8) = Hamiltonian(7,8) - imag*lambda_SOC/2.
    Hamiltonian(10,11) = Hamiltonian(10,11) - imag*lambda_SOC/2.
END SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS

SUBROUTINE COMPUTE_TRIGONAL_TERMS(Hamiltonian, TBA_DIM, DELTA_TRI)
    IMPLICIT NONE 
    INTEGER*4, INTENT(IN):: TBA_DIM
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(2*TBA_DIM,2*TBA_DIM)
    REAL*8, INTENT(IN) :: DELTA_TRI

    !Spin up
    !Ti 1 atoms
    Hamiltonian(1,2) = Hamiltonian(1,2) + DELTA_TRI/2.
    Hamiltonian(1,3) = Hamiltonian(1,3) + DELTA_TRI/2.
    Hamiltonian(2,3) = Hamiltonian(2,3) + DELTA_TRI/2.
    !Ti 2 atoms
    Hamiltonian(1 + 3,2 + 3) = Hamiltonian(1 + 3,2 + 3) + DELTA_TRI/2.
    Hamiltonian(1 + 3,3 + 3) = Hamiltonian(1 + 3,3 + 3) + DELTA_TRI/2.
    Hamiltonian(2 + 3,3 + 3) = Hamiltonian(2 + 3,3 + 3) + DELTA_TRI/2.
    
    !Spin down
    !Ti 1 atoms
    Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) + DELTA_TRI/2.
    !Ti 2 atoms
    Hamiltonian(1 + 3 + TBA_DIM,2 + 3 + TBA_DIM) = Hamiltonian(1 + 3 + TBA_DIM,2 + 3 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(1 + 3 + TBA_DIM,3 + 3 + TBA_DIM) = Hamiltonian(1 + 3 + TBA_DIM,3 + 3 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(2 + 3 + TBA_DIM,3 + 3 + TBA_DIM) = Hamiltonian(2 + 3 + TBA_DIM,3 + 3 + TBA_DIM) + DELTA_TRI/2.
END SUBROUTINE COMPUTE_TRIGONAL_TERMS

SUBROUTINE DIAGONALIZE_HERMITIAN(Hamiltonian, Eigenvalues, dim)
    IMPLICIT NONE 
    INTEGER*4 :: dim
    COMPLEX*16 :: Hamiltonian(dim, dim)
    REAL*8 :: Eigenvalues(dim)
    COMPLEX*16, ALLOCATABLE :: WORK(:)
    REAL*8, ALLOCATABLE :: RWORK(:)
    INTEGER*4 :: LWORK 
    INTEGER*4 :: INFO
    LWORK = 10*dim
    ALLOCATE(WORK(LWORK))
    ALLOCATE(RWORK(3*dim - 2))
    CALL ZHEEV('V', 'U', dim, Hamiltonian, dim, Eigenvalues, WORK, LWORK, RWORK, INFO)
        !WRITE(*, fmt = "(15F10.4)") Eigenvalues
    IF (INFO .ne. 0) THEN 
        PRINT*, 'ZHEEV INFO ', INFO
        STOP
    END IF 

    DEALLOCATE(WORK)
    DEALLOCATE(RWORK)

END SUBROUTINE DIAGONALIZE_HERMITIAN


SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, dim)
    IMPLICIT NONE 
    INTEGER*4 :: dim
    COMPLEX*16 :: Hamiltonian(dim,dim)
    INTEGER*4 :: i,j
    DO i = 1, dim
        DO j = i + 1, dim
            Hamiltonian(j,i) = CONJG(Hamiltonian(i,j))
        END DO
    END DO
END SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS


COMPLEX*16 FUNCTION epsilon_yz(t_D, t_I, kx, ky)
    REAL*8, INTENT(IN) :: t_D, t_I, kx, ky
    epsilon_yz = -t_D*(1. + EXP(DCMPLX(0. , SQRT(3.)/2.*kx - 3./2.*ky))) &
                - t_I*EXP(DCMPLX(0. , -SQRT(3.)/2.*kx - 3./2.*ky ))
    RETURN 
END FUNCTION epsilon_yz 

COMPLEX*16 FUNCTION epsilon_zx(t_D, t_I, kx, ky)
    REAL*8, INTENT(IN) :: t_D, t_I, kx, ky
    epsilon_zx = -t_D*(1. + EXP(DCMPLX(0. , -SQRT(3.)/2.*kx - 3./2.*ky))) &
                - t_I*EXP(DCMPLX(0. , SQRT(3.)/2.*kx - 3./2.*ky ))
    RETURN
END FUNCTION epsilon_zx

COMPLEX*16 FUNCTION epsilon_xy(t_D, t_I, kx, ky)
    REAL*8, INTENT(IN) :: t_D, t_I, kx, ky
    epsilon_xy = -2.*t_D*DCOS( SQRT(3.)/2.*kx )*EXP(DCMPLX(0. -3./2.*ky)) - t_I
    RETURN
END FUNCTION epsilon_xy


REAL*8 FUNCTION meV_to_au(x)
    IMPLICIT NONE
    REAL*8 :: x
    meV_to_au = x / 27211.
    RETURN 
END FUNCTION meV_to_au

REAL*8 FUNCTION nm_to_au(x)
    IMPLICIT NONE
    REAL*8 :: x 
    nm_to_au = x/0.05292
    RETURN
END FUNCTION nm_to_au

REAL*8 FUNCTION au_to_meV(x)
    IMPLICIT NONE
    REAL*8 :: x
    au_to_meV = x * 27211
    RETURN 
END FUNCTION au_to_meV

REAL*8 FUNCTION au_to_nm(x)
    IMPLICIT NONE
    REAL*8 :: x 
    au_to_nm = x*0.05292
    RETURN
END FUNCTION au_to_nm
