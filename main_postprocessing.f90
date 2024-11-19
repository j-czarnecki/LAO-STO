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

    zeta = 1e-3


    CALL GET_POSTPROCESSING_INPUT('./postprocessing_input.nml')

    CALL CALCULATE_DISPERSION("Energies.dat", 3)
    CALL CALCULATE_DOS(-2.1e3*meV2au, 2.1e3*meV2au, 1.0e3*meV2au, zeta, "OutputData/DOS.dat")

    IF (enable_chern_number_calc) THEN
        CALL CALCULATE_CHERN_PARAMS(Nk_points_chern_number, Nk_points_chern_number, path_to_run_dir_chern_number)
    END IF

    IF (enable_sc_gap_calc) THEN
        CALL CALCULATE_SUPERCONDUCTING_GAP(TRIM(path_to_run_dir_sc_gap), dE_sc_gap, Nk_points_sc_gap)
    END IF



END PROGRAM main_postprocessing