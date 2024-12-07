PROGRAM main_postprocessing
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_compute_hamiltonians
    USE mod_postprocessing
    IMPLICIT NONE


    CALL GET_POSTPROCESSING_INPUT('./postprocessing_input.nml')


    IF (enable_dispersion_relation_calc) THEN
        CALL CALCULATE_DISPERSION(path_to_run_dir_dispersion_relation, Nk_points_dispersion_relation)
    END IF

    IF (enable_dos_calc) THEN
        CALL CALCULATE_DOS(E_DOS_min, E_DOS_max, dE0, zeta_DOS, Nk_points_dos, path_to_run_dir_dos)
    END IF

    IF (enable_chern_number_calc) THEN
        CALL CALCULATE_CHERN_PARAMS(Nk_points_chern_number, Nk_points_chern_number, path_to_run_dir_chern_number)
    END IF

    IF (enable_sc_gap_calc) THEN
        CALL CALCULATE_SUPERCONDUCTING_GAP(TRIM(path_to_run_dir_sc_gap), dE_sc_gap, Nk_points_sc_gap)
    END IF



END PROGRAM main_postprocessing