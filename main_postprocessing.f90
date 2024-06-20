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
    INTEGER*4 :: Nk, HamDim

    zeta = 1e-3

    Nk = 5000
    HamDim = DIM


    !CALL CALCULATE_DISPERSION("Energies.dat")
    CALL CALCULATE_DOS(-1.1e3*meV2au, -0.9e3*meV2au, 1.*meV2au, zeta, "OutputData/DOS.dat")
    !CALL CALCULATE_CHERN_PARAMS(Nk, Nk, HamDim)


END PROGRAM main_postprocessing