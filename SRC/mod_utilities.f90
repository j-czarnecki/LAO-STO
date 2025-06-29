MODULE mod_utilities
USE mod_parameters
USE mod_reader
IMPLICIT NONE
CONTAINS

SUBROUTINE DIAGONALIZE_HERMITIAN(Hamiltonian, Eigenvalues, N)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: N
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(N, N)
  REAL*8, INTENT(OUT) :: Eigenvalues(N)
  COMPLEX*16, ALLOCATABLE :: WORK(:)
  REAL*8, ALLOCATABLE :: RWORK(:)
  INTEGER*4 :: LWORK
  INTEGER*4 :: INFO
  LWORK = 10 * DIM
  ALLOCATE (WORK(LWORK))
  ALLOCATE (RWORK(3 * N - 2))
  CALL ZHEEV('V', 'U', N, Hamiltonian, N, Eigenvalues, WORK, LWORK, RWORK, INFO)
  !WRITE(*, fmt = "(15F10.4)") Eigenvalues
  IF (INFO .ne. 0) THEN
    PRINT *, 'ZHEEV INFO ', INFO
    STOP
  END IF

  DEALLOCATE (WORK)
  DEALLOCATE (RWORK)

END SUBROUTINE DIAGONALIZE_HERMITIAN

SUBROUTINE DIAGONALIZE_GENERALIZED(Hamiltonian, Eigenvalues, U_transformation, N)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: N
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(N, N)
  COMPLEX*16, INTENT(INOUT) :: U_transformation(N, N)
  REAL*8, INTENT(OUT) :: Eigenvalues(N)
  COMPLEX*16, ALLOCATABLE :: W(:)
  COMPLEX*16, ALLOCATABLE :: VL(:, :)
  COMPLEX*16, ALLOCATABLE :: WORK(:)
  REAL*8, ALLOCATABLE :: RWORK(:)
  INTEGER*4 :: LWORK
  INTEGER*4 :: INFO, i

  LWORK = 10 * N
  INFO = 0

  ALLOCATE (W(N))
  ALLOCATE (VL(N, N))
  ALLOCATE (WORK(LWORK))
  ALLOCATE (RWORK(2 * N))

  W(:) = 0.
  VL(:, :) = DCMPLX(0., 0.)
  WORK(:) = 0.
  RWORK(:) = 0.

  CALL ZGEEV('N', 'V', N, Hamiltonian, N, W, VL, N, U_transformation, N,&
  & WORK, LWORK, RWORK, INFO)

  IF (INFO .ne. 0) THEN
    PRINT *, 'ZGEEV INFO ', INFO
    STOP
  END IF

  !Removing phase ambiguity
  DO i = 1, N
    U_transformation(:, i) = U_transformation(:, i) * CONJG(U_transformation(1, i)) / ABS(U_transformation(1, i))
  END DO
  Eigenvalues(:) = REAL(W(:))

  DEALLOCATE (W)
  DEALLOCATE (VL)
  DEALLOCATE (WORK)
  DEALLOCATE (RWORK)

END SUBROUTINE DIAGONALIZE_GENERALIZED

SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS(Hamiltonian, N)
  IMPLICIT NONE
  INTEGER*4, INTENT(IN) :: N
  COMPLEX*16, INTENT(INOUT) :: Hamiltonian(N, N)
  INTEGER*4 :: i, j
  DO i = 1, N
    DO j = i + 1, N
      Hamiltonian(j, i) = CONJG(Hamiltonian(i, j))
    END DO
  END DO
END SUBROUTINE COMPUTE_CONJUGATE_ELEMENTS

!dir$ attributes forceinline :: epsilon_yz
PURE COMPLEX * 16 FUNCTION epsilon_yz(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  epsilon_yz = -t_D * (EXP(imag * ky) + EXP(imag * (SQRT(3.) / 2.*kx - 1./2.*ky))) &
               - t_I * EXP(-imag * (SQRT(3.) / 2.*kx + 1./2.*ky))
  ! epsilon_yz = -t_D * EXP(-imag * ky) * (1.0 + EXP(-imag * (SQRT(3.) / 2.*kx - 3./2.*ky))) &
  !              - t_I * EXP(-imag * (-SQRT(3.) / 2.*kx - 1./2.*ky))
  RETURN
END FUNCTION epsilon_yz

!dir$ attributes forceinline :: epsilon_zx
PURE COMPLEX * 16 FUNCTION epsilon_zx(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  epsilon_zx = -t_D * (EXP(imag * ky) + EXP(-imag * (SQRT(3.) / 2.*kx + 1./2.*ky))) &
               - t_I * EXP(imag * (SQRT(3.) / 2.*kx - 1./2.*ky))
  ! epsilon_zx = -t_D * EXP(-imag * ky) * (1.+EXP(-imag * (-SQRT(3.) / 2.*kx - 3./2.*ky))) &
  !              - t_I * EXP(-imag * (SQRT(3.) / 2.*kx - 1./2.*ky))
  RETURN
END FUNCTION epsilon_zx

!dir$ attributes forceinline :: epsilon_xy
PURE COMPLEX * 16 FUNCTION epsilon_xy(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  epsilon_xy = -2.*t_D * COS(SQRT(3.) / 2.*kx) * EXP(-imag * 1./2.*ky) - t_I * EXP(imag * ky)
  ! epsilon_xy = -2.*t_D * COS(SQRT(3.) / 2.*kx) * EXP(imag * 1./2.*ky) - t_I * EXP(-imag * ky)
  RETURN
END FUNCTION epsilon_xy

!dir$ attributes forceinline :: rashba_yz_xz
PURE COMPLEX * 16 FUNCTION rashba_yz_zx(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  rashba_yz_zx = 2 * imag * t_Rashba * SIN(-SQRT(3.) / 2.*kx) * EXP(-imag * 1./2.*ky)
  RETURN
END FUNCTION rashba_yz_zx

!dir$ attributes forceinline :: rashba_yz_xy
PURE COMPLEX * 16 FUNCTION rashba_yz_xy(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  rashba_yz_xy = -t_Rashba * EXP(imag * ky) * (1.-EXP(imag * (-SQRT(3.) / 2.*kx - 3./2.*ky)))
  RETURN
END FUNCTION rashba_yz_xy

!dir$ attributes forceinline :: rashba_zx_xy
PURE COMPLEX * 16 FUNCTION rashba_zx_xy(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  rashba_zx_xy = -t_Rashba * EXP(imag * ky) * (1.-EXP(imag * (SQRT(3.) / 2.*kx - 3./2.*ky)))
  RETURN
END FUNCTION rashba_zx_xy

!dir$ attributes forceinline :: pairing_1
PURE COMPLEX * 16 FUNCTION pairing_1(ky)
  REAL*8, INTENT(IN) :: ky
  pairing_1 = EXP(imag * ky)
  RETURN
END FUNCTION pairing_1

!dir$ attributes forceinline :: pairing_2
PURE COMPLEX * 16 FUNCTION pairing_2(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_2 = EXP(-imag * (SQRT(3.) / 2.*kx + ky / 2.))
  RETURN
END FUNCTION pairing_2

!dir$ attributes forceinline :: pairing_3
PURE COMPLEX * 16 FUNCTION pairing_3(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_3 = EXP(-imag * (-SQRT(3.) / 2.*kx + ky / 2.))
  RETURN
END FUNCTION pairing_3

!dir$ attributes forceinline :: pairing_nnn_1
PURE COMPLEX * 16 FUNCTION pairing_nnn_1(kx)
  REAL*8, INTENT(IN) :: kx
  pairing_nnn_1 = EXP(-imag * (SQRT(3.) * kx))
  RETURN
END FUNCTION pairing_nnn_1

!dir$ attributes forceinline :: pairing_nnn_2
PURE COMPLEX * 16 FUNCTION pairing_nnn_2(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_nnn_2 = EXP(-imag * (SQRT(3.) / 2.*kx + 3./2.*ky))
  RETURN
END FUNCTION pairing_nnn_2

!dir$ attributes forceinline :: pairing_nnn_3
PURE COMPLEX * 16 FUNCTION pairing_nnn_3(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_nnn_3 = EXP(-imag * (-SQRT(3.) / 2.*kx + 3./2.*ky))
  RETURN
END FUNCTION pairing_nnn_3

!dir$ attributes forceinline :: pairing_nnn_4
PURE COMPLEX * 16 FUNCTION pairing_nnn_4(kx)
  REAL*8, INTENT(IN) :: kx
  pairing_nnn_4 = EXP(-imag * (-SQRT(3.) * kx))
  RETURN
END FUNCTION pairing_nnn_4

!dir$ attributes forceinline :: pairing_nnn_5
PURE COMPLEX * 16 FUNCTION pairing_nnn_5(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_nnn_5 = EXP(-imag * (-SQRT(3.) / 2.*kx - 3./2.*ky))
  RETURN
END FUNCTION pairing_nnn_5

!dir$ attributes forceinline :: pairing_nnn_6
PURE COMPLEX * 16 FUNCTION pairing_nnn_6(kx, ky)
  REAL*8, INTENT(IN) :: kx, ky
  pairing_nnn_6 = EXP(-imag * (SQRT(3.) / 2.*kx - 3./2.*ky))
  RETURN
END FUNCTION pairing_nnn_6

!dir$ attributes forceinline :: fd_distribution
RECURSIVE PURE REAL * 8 FUNCTION fd_distribution(E, E_Fermi, T)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: E, E_Fermi, T
  fd_distribution = 0.
  IF (T .ne. 0.0d0) THEN
    fd_distribution = 1./(EXP((E - E_Fermi) / (k_B * T)) + 1.)
  ELSE
    IF (E > E_Fermi) THEN
      fd_distribution = 0.
    ELSE IF (E == E_Fermi) THEN
      fd_distribution = 0.5
    ELSE IF (E < E_Fermi) THEN
      fd_distribution = 1.
    END IF
  END IF
  RETURN
END FUNCTION fd_distribution

!dir$ attributes forceinline :: fd_distribution_derivative_E
RECURSIVE PURE REAL * 8 FUNCTION fd_distribution_derivative_E(E, E_Fermi, T)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: E, E_Fermi, T
  REAL*8 :: x
  x = EXP((E - E_Fermi) / (k_B * T))
  fd_distribution_derivative_E = -x / ((k_B * T) * (x + 1.0d0)**2)
  RETURN
END FUNCTION fd_distribution_derivative_E

!dir$ attributes forceinline :: dirac_delta
RECURSIVE PURE REAL * 8 FUNCTION dirac_delta(E, omega, zeta)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: E, omega, zeta
  dirac_delta = zeta / (PI * ((E - omega)**2 + zeta**2))
  RETURN
END FUNCTION dirac_delta

!dir$ attributes forceinline :: r_max_phi
RECURSIVE PURE REAL * 8 FUNCTION r_max_phi(phi)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: phi
  r_max_phi = R_K_MAX * SQRT(3.0d0) / (2.0d0 * COS(phi - PI / 6.0d0))
  RETURN
END FUNCTION r_max_phi

END MODULE mod_utilities
