#include "macros_def.f90"
PROGRAM main_postprocessing
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_utilities
    USE mod_writers
    USE mod_reader
    USE mod_compute_hamiltonians
    USE mod_postprocessing
    USE mod_logger
    IMPLICIT NONE

    CALL INIT_LOGGER("postprocessing")

    WRITE (log_string,*) "Reading postprocessing input"
    LOG_INFO(log_string)

    CALL GET_POSTPROCESSING_INPUT('./postprocessing_input.nml')


    IF (enable_dispersion_relation_calc) THEN
        WRITE (log_string,*) "Calculating dispersion relation"
        LOG_INFO(log_string)
        CALL CALCULATE_DISPERSION(path_to_run_dir_dispersion_relation, Nk_points_dispersion_relation, include_sc_in_dispersion)
    END IF

    IF (enable_dos_calc) THEN
        WRITE (log_string,*) "Calculating DOS"
        LOG_INFO(log_string)
        CALL CALCULATE_DOS(E_DOS_min, E_DOS_max, dE0, zeta_DOS, include_sc_in_dos, Nk_points_dos, path_to_run_dir_dos)
    END IF

    IF (enable_chern_number_calc) THEN
        WRITE (log_string,*) "Calculating Chern number"
        LOG_INFO(log_string)
        CALL CALCULATE_CHERN_PARAMS(Nk_points_chern_number, Nk_points_chern_number, path_to_run_dir_chern_number)
    END IF

    IF (enable_sc_gap_calc) THEN
        WRITE (log_string,*) "Calculating superconducting gap"
        LOG_INFO(log_string)
        CALL CALCULATE_SUPERCONDUCTING_GAP(TRIM(path_to_run_dir_sc_gap), dE_sc_gap, Nk_points_sc_gap)
    END IF

    CALL CLOSE_LOGGER()

    IF (ALLOCATED(V_layer)) DEALLOCATE(V_layer)

END PROGRAM main_postprocessing