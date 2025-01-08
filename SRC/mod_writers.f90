MODULE mod_writers
USE mod_parameters
USE mod_reader
IMPLICIT NONE
CONTAINS

SUBROUTINE PRINT_HAMILTONIAN(Hamiltonian)
    COMPLEX*16, INTENT(IN) :: Hamiltonian(DIM,DIM)
    CHARACTER(LEN=20) :: output_format
    INTEGER*4 :: i
    
    WRITE(output_format, '(A, I0, A)') '(',DIM,'E15.5)'
    output_format = TRIM(output_format)

    OPEN(unit = 9, FILE= "./OutputData/H_real.dat", FORM = "FORMATTED", ACTION = "WRITE")
    OPEN(unit = 10, FILE= "./OutputData/H_imag.dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO i = 1, DIM
        WRITE(9, output_format) REAL(Hamiltonian(i,:))
        WRITE(10, output_format) AIMAG(Hamiltonian(i,:))
    END DO
    CLOSE(9)
    CLOSE(10)


END SUBROUTINE PRINT_HAMILTONIAN

SUBROUTINE PRINT_ENERGIES(Energies, k1_steps, k2_steps, dk1, dk2, filename, N)
    INTEGER*4, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: Energies(0:k1_steps, 0:k2_steps, N)
    REAL*8, INTENT(IN) :: dk1, dk2
    INTEGER*4, INTENT(IN) :: k1_steps, k2_steps
    REAL*8 :: k1, k2, kx, ky

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format
    INTEGER*4 :: i,j,l

    output_format = '(I5, 3E15.5)'

    OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
    DO l = 1, N
        DO i = 0,k1_steps
            DO j = 0, k2_steps
                k1 = i*dk1
                k2 = j*dk2

                kx = 2.*PI/(SQRT(3.0d0)) * k1
                ky = -2.*PI/3. * k1 + 4.*PI/3. * k2
                WRITE(9, output_format) l, k1, k2, Energies(i, j, l)/meV2au
            END DO
            WRITE(9,*)
            WRITE(9,*)
        END DO
        WRITE(9,*)
        WRITE(9,*)
    END DO
    CLOSE(9)
END SUBROUTINE

SUBROUTINE PRINT_GAMMA(Gamma_SC, filename)
    COMPLEX*16, INTENT(IN) :: Gamma_SC(ORBITALS,N_ALL_NEIGHBOURS, 2,LAYER_COUPLINGS, SUBBANDS)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format

    INTEGER*4 :: orb, j,spin, lat, band
    output_format = '(5I5, 2E15.5)'

    !Printing SC gammas in [meV]
    OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,*) "#band spin neighbour lattice orbital Re(Gamma) Im(Gamma)"
    DO band = 1, SUBBANDS
        DO spin = 1, 2
            !Do this, because we store LAYER_COUPLINGS for nearest neighbours, but only SUBLATTICES for next-to-nearest
            DO j = 1, N_NEIGHBOURS
                DO lat = 1, LAYER_COUPLINGS
                    DO orb = 1, ORBITALS
                        WRITE(9, output_format) band, spin, j, lat, orb, REAL(Gamma_SC(orb,j, spin,lat, band))/meV2au, AIMAG(Gamma_SC(orb,j, spin,lat, band))/meV2au
                    END DO
                END DO
                WRITE(9,*)
                WRITE(9,*)
            END DO
            DO j = N_NEIGHBOURS + 1, N_ALL_NEIGHBOURS
                DO lat = 1, SUBLATTICES
                    DO orb = 1, ORBITALS
                        WRITE(9, output_format) band, spin, j, lat, orb, REAL(Gamma_SC(orb,j, spin,lat, band))/meV2au, AIMAG(Gamma_SC(orb,j, spin,lat, band))/meV2au
                    END DO
                END DO
                WRITE(9,*)
                WRITE(9,*)
            END DO

        END DO
    END DO
    CLOSE(9)
END SUBROUTINE PRINT_GAMMA

SUBROUTINE PRINT_CHARGE(Charge_dens, filename)
    REAL*8, INTENT(IN) :: Charge_dens(DIM_POSITIVE_K, SUBBANDS)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=20) :: output_format
    INTEGER*4 :: spin, lat, orb, n, band

    output_format = '(4I5, 1E15.5)'
    OPEN(unit = 9, FILE= "./OutputData/"//filename//".dat", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(9,*) "#band spin lattice orbital Charge"
    DO band = 1, SUBBANDS
        n = 1
        DO spin = 1, 2
            DO lat = 1, SUBLATTICES
                DO orb = 1, ORBITALS
                    WRITE(9, output_format) band, spin, lat, orb, Charge_dens(n, band)
                    n = n + 1
                END DO
            END DO
        END DO
    END DO
    CLOSE(9)
END SUBROUTINE PRINT_CHARGE

END MODULE mod_writers