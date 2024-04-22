PROGRAM main_postprocessing
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_compute_hamiltonians
    USE mod_postprocessing
    IMPLICIT NONE

    REAL*8 :: zeta

    zeta = 1e-5

    CALL CALCULATE_DISPERSION("Energies.dat")
    CALL CALCULATE_DOS(-1.2e3*meV2au, -0.8e3*meV2au, 0.1*meV2au, zeta, "OutputData/DOS.dat")

    
    !CALL CALCULATE_CHERN()

END PROGRAM main_postprocessing