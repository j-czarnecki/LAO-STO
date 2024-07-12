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

    zeta = 1e-6

    Nk = 5000
    HamDim = DIM


    CALL CALCULATE_DISPERSION("Energies.dat", 3)
    !CALL CALCULATE_DOS(-1.1e3*meV2au, -0.9e3*meV2au, 0.01*meV2au, zeta, "OutputData/DOS.dat")
    !CALL CALCULATE_CHERN_PARAMS(Nk, Nk, HamDim)
    ! CALL CALCULATE_SUPERCONDUCTING_GAP("/home/jczarnecki/LAO-STO-results/LAO-STO/RUN_E_Fermi_-1000.0_J_SC_150.0/", "_NN_J150", 0.6d0, 2e-3*meV2au, 10000)
    ! PRINT*, 'Done'
    ! CALL CALCULATE_SUPERCONDUCTING_GAP("/home/jczarnecki/LAO-STO-results/LAO-STO/RUN_E_Fermi_-1000.0/", "_NN_J175", 0.6d0, 2e-3*meV2au, 10000)
    ! PRINT*, 'Done'
    ! CALL CALCULATE_SUPERCONDUCTING_GAP("/home/jczarnecki/LAO-STO-results/LAO-STO-nnn/RUN_E_Fermi_-1000.0_J_SC_NNN_125.0/", "_NNN_J150", 0.6d0, 2e-3*meV2au, 10000)
    ! PRINT*, 'Done'
    ! CALL CALCULATE_SUPERCONDUCTING_GAP("/home/jczarnecki/LAO-STO-results/LAO-STO-nnn/RUN_E_Fermi_-1000.0/", "_NNN_J175", 0.6d0, 2e-3*meV2au, 10000)
    ! PRINT*, 'Done'

    !CALL CALCULATE_SUPERCONDUCTING_GAP("./", "", 0.2d0, 5e-4*meV2au, 12000)

END PROGRAM main_postprocessing