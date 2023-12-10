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










END MODULE mod_writers