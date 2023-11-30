MODULE mod_hamiltonians
USE mod_utilities
USE mod_parameters
IMPLICIT NONE 
CONTAINS


SUBROUTINE COMPUTE_TBA_TERM(Hamiltonian, kx, ky)
    IMPLICIT NONE 
    REAL*8 :: kx, ky
    COMPLEX*16 :: Hamiltonian(DIM,DIM) !Twice as big because of spin
    !Only specifying upper triangle of matrix, since Hamiltonian is hermitian
    !spin up
    Hamiltonian(1,4) = epsilon_yz(kx, ky) 
    Hamiltonian(2,5) = epsilon_zx(kx, ky) 
    Hamiltonian(3,6) = epsilon_xy(kx, ky) 
    !spin down
    Hamiltonian(TBA_DIM + 1,TBA_DIM + 4) = epsilon_yz(kx, ky)
    Hamiltonian(TBA_DIM + 2,TBA_DIM + 5) = epsilon_zx(kx, ky)
    Hamiltonian(TBA_DIM + 3,TBA_DIM + 6) = epsilon_xy(kx, ky)
END SUBROUTINE COMPUTE_TBA_TERM

SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS(Hamiltonian)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)

    !Ti1 atoms
    Hamiltonian(1,2) = Hamiltonian(1,2) + imag*lambda_SOC/2.
    Hamiltonian(7,8) = Hamiltonian(7,8) - imag*lambda_SOC/2.
    Hamiltonian(1,9) = Hamiltonian(1,9) - lambda_SOC/2.
    Hamiltonian(2,9) = Hamiltonian(2,9) + imag*lambda_SOC/2.
    Hamiltonian(3,7) = Hamiltonian(3,7) + lambda_SOC/2.
    Hamiltonian(3,8) = Hamiltonian(3,8) - imag*lambda_SOC/2.

    !Ti2 atoms
    Hamiltonian(4,5) = Hamiltonian(4,5) + imag*lambda_SOC/2.
    Hamiltonian(10,11) = Hamiltonian(10,11) - imag*lambda_SOC/2.
    Hamiltonian(4,12) = Hamiltonian(4,12) - lambda_SOC/2.
    Hamiltonian(5,12) = Hamiltonian(5,12) + imag*lambda_SOC/2.
    Hamiltonian(6,10) = Hamiltonian(6,10) + lambda_SOC/2.
    Hamiltonian(6,11) = Hamiltonian(6,11) - imag*lambda_SOC/2.

END SUBROUTINE COMPUTE_ATOMIC_SOC_TERMS

SUBROUTINE COMPUTE_TRIGONAL_TERMS(Hamiltonian)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(2*TBA_DIM,2*TBA_DIM)

    !Spin up
    !Ti 1 atoms
    Hamiltonian(1,2) = Hamiltonian(1,2) + DELTA_TRI/2.
    Hamiltonian(1,3) = Hamiltonian(1,3) + DELTA_TRI/2.
    Hamiltonian(2,3) = Hamiltonian(2,3) + DELTA_TRI/2.
    !Ti 2 atoms
    Hamiltonian(ORBITALS + 1, ORBITALS + 2) = Hamiltonian(ORBITALS + 1, ORBITALS + 2) + DELTA_TRI/2.
    Hamiltonian(ORBITALS + 1, ORBITALS + 3) = Hamiltonian(ORBITALS + 1, ORBITALS + 3) + DELTA_TRI/2.
    Hamiltonian(ORBITALS + 2, ORBITALS + 3) = Hamiltonian(ORBITALS + 2, ORBITALS + 3) + DELTA_TRI/2.
    
    !Spin down
    !Ti 1 atoms
    Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) + DELTA_TRI/2.
    !Ti 2 atoms
    Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) + DELTA_TRI/2.
    Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) + DELTA_TRI/2.
END SUBROUTINE COMPUTE_TRIGONAL_TERMS


SUBROUTINE COMPUTE_ELECTRIC_FIELD(Hamiltonian)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)
    INTEGER*4 :: i

    DO i = 1, 3
        !Ti1 atoms
        Hamiltonian(i,i) = Hamiltonian(i,i) + v/2.
        Hamiltonian(i + TBA_DIM, i + TBA_DIM) = Hamiltonian(i + TBA_DIM, i + TBA_DIM) + v/2.
        !Ti2 atoms
        Hamiltonian(i + ORBITALS, i + ORBITALS) = Hamiltonian(i + ORBITALS, i + ORBITALS) - v/2.
        Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) = Hamiltonian(i + TBA_DIM + ORBITALS, i + TBA_DIM + ORBITALS) - v/2
    END DO
END SUBROUTINE COMPUTE_ELECTRIC_FIELD

SUBROUTINE COMPUTE_TI1_TI2(Hamiltonian,kx, ky)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)
    REAL*8 :: kx, ky
    !Spin-up part
    Hamiltonian(1,2 + ORBITALS) = Hamiltonian(1, 2 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 2.*imag*EXP(imag*3./2.*ky)*SIN(SQRT(3.)/2.*kx))
    Hamiltonian(1,3 + ORBITALS) = Hamiltonian(1, 3 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(1 - EXP(imag/2.*(SQRT(3.)*kx + 3.*ky)))
    Hamiltonian(2, 1 + ORBITALS) = Hamiltonian(2, 1 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(2*imag*EXP(imag*3./2.*ky)*SIN(SQRT(3.)/2.*kx))
    Hamiltonian(2, 3 + ORBITALS) = Hamiltonian(2, 3 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(1 - EXP(-imag/2.*(SQRT(3.)*kx - 3.*ky)))
    Hamiltonian(3, 1 + ORBITALS) = Hamiltonian(3, 1 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 1 + EXP(imag/2.*(SQRT(3.)*kx + 3.*ky)))
    Hamiltonian(3, 2 + ORBITALS) = Hamiltonian(3, 2 + ORBITALS) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 1 + EXP(-imag/2.*(SQRT(3.)*kx - 3.*ky)))

    !Spin-down part
    Hamiltonian(1 + TBA_DIM,2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 2.*imag*EXP(imag*3./2.*ky)*SIN(SQRT(3.)/2.*kx))
    Hamiltonian(1 + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(1 - EXP(imag/2.*(SQRT(3.)*kx + 3.*ky)))
    Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(2*imag*EXP(imag*3./2.*ky)*SIN(SQRT(3.)/2.*kx))
    Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + TBA_DIM, 3 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(1 - EXP(-imag/2.*(SQRT(3.)*kx - 3.*ky)))
    Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 1 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 1 + EXP(imag/2.*(SQRT(3.)*kx + 3.*ky)))
    Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) = Hamiltonian(3 + TBA_DIM, 2 + ORBITALS + TBA_DIM) + eta_p*V_pdp*SQRT(2.)**(7./4.)/SQRT(15.)*(- 1 + EXP(-imag/2.*(SQRT(3.)*kx - 3.*ky)))

END SUBROUTINE COMPUTE_TI1_TI2

SUBROUTINE COMPUTE_H_PI(Hamiltonian, kx, ky)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)
    REAL*8 :: kx, ky
    REAL*8 :: k1, k2, k3

    k1 = -SQRT(3.)/2.*kx + 3./2.*ky 
    k2 = -SQRT(3.)/2.*kx - 3./2.*ky 
    k3 = SQRT(3.)*kx

    !Spin up, Ti1
    Hamiltonian(1,2) = Hamiltonian(1,2) - eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + SIN(k2) + 2*SIN(k3))
    Hamiltonian(1,3) = Hamiltonian(1,3) + eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + 2*SIN(k2) + SIN(k3))
    Hamiltonian(2,3) = Hamiltonian(2,3) - eta_p*2*imag*V_pdp/SQRT(15.)*(2*SIN(k1) + SIN(k2) + SIN(k3))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS,2 + ORBITALS) = Hamiltonian(1 + ORBITALS,2 + ORBITALS) - eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + SIN(k2) + 2*SIN(k3))
    Hamiltonian(1 + ORBITALS,3 + ORBITALS) = Hamiltonian(1 + ORBITALS,3 + ORBITALS) + eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + 2*SIN(k2) + SIN(k3))
    Hamiltonian(2 + ORBITALS,3 + ORBITALS) = Hamiltonian(2 + ORBITALS,3 + ORBITALS) - eta_p*2*imag*V_pdp/SQRT(15.)*(2*SIN(k1) + SIN(k2) + SIN(k3))

    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) - eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + SIN(k2) + 2*SIN(k3))
    Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) + eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + 2*SIN(k2) + SIN(k3))
    Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) - eta_p*2*imag*V_pdp/SQRT(15.)*(2*SIN(k1) + SIN(k2) + SIN(k3))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) - eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + SIN(k2) + 2*SIN(k3))
    Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) + eta_p*2*imag*V_pdp/SQRT(15.)*(SIN(k1) + 2*SIN(k2) + SIN(k3))
    Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) - eta_p*2*imag*V_pdp/SQRT(15.)*(2*SIN(k1) + SIN(k2) + SIN(k3))

END SUBROUTINE COMPUTE_H_PI

SUBROUTINE COMPUTE_H_SIGMA(Hamiltonian, kx, ky)
    IMPLICIT NONE 
    COMPLEX*16, INTENT(INOUT) :: Hamiltonian(DIM,DIM)
    REAL*8 :: kx, ky
    REAL*8 :: k1, k2, k3
    k1 = -SQRT(3.)/2.*kx + 3./2.*ky 
    k2 = -SQRT(3.)/2.*kx - 3./2.*ky 
    k3 = SQRT(3.)*kx

    !Spin up, Ti1
    Hamiltonian(1,2) = Hamiltonian(1,2) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(1,3) = Hamiltonian(1,3) - eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(2,3) = Hamiltonian(2,3) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k2) + SIN(k3))
    !Spin up, Ti2
    Hamiltonian(1 + ORBITALS,2 + ORBITALS) = Hamiltonian(1 + ORBITALS,2 + ORBITALS) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(1 + ORBITALS,3 + ORBITALS) = Hamiltonian(1 + ORBITALS,3 + ORBITALS) - eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(2 + ORBITALS,3 + ORBITALS) = Hamiltonian(2 + ORBITALS,3 + ORBITALS) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k2) + SIN(k3))
    !Spin down, Ti1
    Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,2 + TBA_DIM) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(1 + TBA_DIM,3 + TBA_DIM) - eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) = Hamiltonian(2 + TBA_DIM,3 + TBA_DIM) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k2) + SIN(k3))
    !Spin down, Ti2
    Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,2 + ORBITALS + TBA_DIM) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(1 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) - eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k1) + SIN(k2))
    Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) = Hamiltonian(2 + ORBITALS + TBA_DIM,3 + ORBITALS + TBA_DIM) + eta_p*2*imag*SQRT(3.)*V_pds/SQRT(15.)*(SIN(k2) + SIN(k3))

END SUBROUTINE COMPUTE_H_SIGMA



END MODULE mod_hamiltonians