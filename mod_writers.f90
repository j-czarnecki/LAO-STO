MODULE mod_writers
USE mod_parameters
IMPLICIT NONE
CONTAINS

SUBROUTINE PRINT_HAMILTONIAN(Hamiltonian, filename)
    REAL*8, INTENT(IN) :: Hamiltonian(DIM,DIM)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format
    INTEGER*4 :: i
    
    WRITE(output_format, '(A, I0, A)') '(',DIM,'E15.5)'
    output_format = TRIM(output_format)

    IF(LEN(filename) .NE. 0) THEN
        OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
        DO i = 1, DIM
            WRITE(9, output_format) Hamiltonian(i,:)
        END DO
        CLOSE(9)
    ELSE
        DO i = 1, DIM
            WRITE(*, output_format) Hamiltonian(i,:)
        END DO        
    END IF

END SUBROUTINE PRINT_HAMILTONIAN

SUBROUTINE PRINT_ENERGIES(Energies, k1_steps, k2_steps, dk1, dk2, filename)
    REAL*8, INTENT(IN) :: Energies(0:k1_steps, 0:k2_steps, DIM)
    REAL*8, INTENT(IN) :: dk1, dk2
    INTEGER*4, INTENT(IN) :: k1_steps, k2_steps

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format
    INTEGER*4 :: i,j,n

    output_format = '(E15.5)'

    OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO n = 1, DIM
        DO i = 0,k1_steps
            DO j = 0, k2_steps
                WRITE(9,*) i*dk1*SQRT(3.)/2. * A_TILDE, ( -i*dk1/2. + j*dk2 ) * A_TILDE, Energies(i, j, n)/meV2au*1e-3
            END DO
        END DO
        WRITE(9,*)
        WRITE(9,*)
    END DO 
    CLOSE(9)
END SUBROUTINE

SUBROUTINE PRINT_GAMMA(Gamma_SC, filename)
    COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS,N_NEIGHBOURS, 2)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format

    INTEGER*4 :: orb, j,spin
    output_format = '(3I5, 2E15.5)'

    OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO spin =1, 2
        DO j = 1, N_NEIGHBOURS
            DO orb = 1, ORBITALS
                WRITE(9, output_format) spin, j, orb, REAL(Gamma_SC(orb,j, spin)), AIMAG(Gamma_SC(orb,j, spin))
            END DO
            WRITE(9,*)
            WRITE(9,*)
        END DO
    END DO
    CLOSE(9)
END SUBROUTINE PRINT_GAMMA










END MODULE mod_writers