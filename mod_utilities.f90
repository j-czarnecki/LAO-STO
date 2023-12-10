MODULE mod_utilities 
USE mod_parameters
IMPLICIT NONE
CONTAINS

SUBROUTINE DIAGONALIZE_HERMITIAN(Hamiltonian, Eigenvalues)
    IMPLICIT NONE 
    COMPLEX*16 :: Hamiltonian(DIM, DIM)
    REAL*8 :: Eigenvalues(DIM)
    COMPLEX*16, ALLOCATABLE :: WORK(:)
    REAL*8, ALLOCATABLE :: RWORK(:)
    INTEGER*4 :: LWORK 
    INTEGER*4 :: INFO
    LWORK = 10*DIM
    ALLOCATE(WORK(LWORK))
    ALLOCATE(RWORK(3*DIM - 2))
    CALL ZHEEV('V', 'U', DIM, Hamiltonian, DIM, Eigenvalues, WORK, LWORK, RWORK, INFO)
        !WRITE(*, fmt = "(15F10.4)") Eigenvalues
    IF (INFO .ne. 0) THEN 
        PRINT*, 'ZHEEV INFO ', INFO
        STOP
    END IF 

    DEALLOCATE(WORK)
    DEALLOCATE(RWORK)

END SUBROUTINE DIAGONALIZE_HERMITIAN


SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian)
    IMPLICIT NONE 
    COMPLEX*16 :: Hamiltonian(DIM,DIM)
    INTEGER*4 :: i,j
    DO i = 1, DIM
        DO j = i + 1, DIM
            Hamiltonian(j,i) = CONJG(Hamiltonian(i,j))
        END DO
    END DO
END SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS


COMPLEX*16 FUNCTION epsilon_yz(kx, ky)
    REAL*8, INTENT(IN) :: kx, ky
    epsilon_yz = -t_D*(1. + EXP(imag*(SQRT(3.)/2.*kx - 3./2.*ky))) &
                - t_I*EXP(-imag*(SQRT(3.)/2.*kx + 3./2.*ky ))
    RETURN 
END FUNCTION epsilon_yz 

COMPLEX*16 FUNCTION epsilon_zx(kx, ky)
    REAL*8, INTENT(IN) :: kx, ky
    epsilon_zx = -t_D*(1. + EXP(-imag*(SQRT(3.)/2.*kx + 3./2.*ky))) &
                - t_I*EXP(imag*(SQRT(3.)/2.*kx - 3./2.*ky ))
    RETURN
END FUNCTION epsilon_zx

COMPLEX*16 FUNCTION epsilon_xy(kx, ky)
    REAL*8, INTENT(IN) :: kx, ky
    epsilon_xy = -2.*t_D*DCOS( SQRT(3.)/2.*kx )*EXP(-imag*3./2.*ky) - t_I
    RETURN
END FUNCTION epsilon_xy

COMPLEX*16 FUNCTION pairing_1(ky)
    REAL*8, INTENT(IN) :: ky
    pairing_1 = EXP(-imag*ky)
    RETURN
END FUNCTION pairing_1

COMPLEX*16 FUNCTION pairing_2(kx, ky)
    REAL*8, INTENT(IN) :: kx, ky
    pairing_2 = EXP(imag*(SQRT(3.)/2.*kx + ky/2.))
    RETURN
END FUNCTION pairing_2

COMPLEX*16 FUNCTION pairing_3(kx, ky)
    REAL*8, INTENT(IN) :: kx, ky
    pairing_3 = EXP(imag*(-SQRT(3.)/2.*kx + ky/2.))
    RETURN
END FUNCTION pairing_3

REAL*8 FUNCTION meV_to_au(x)
    IMPLICIT NONE
    REAL*8 :: x
    meV_to_au = x / 27211.
    RETURN 
END FUNCTION meV_to_au

REAL*8 FUNCTION fd_distribution(E, E_Fermi, T)
    IMPLICIT NONE 
    REAL*8 E, E_Fermi, T
    fd_distribution = 1./(EXP((E - E_Fermi)/(k_B*T)) + 1.)
    RETURN
END FUNCTION fd_distribution

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


END MODULE mod_utilities