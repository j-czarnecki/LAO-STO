from RunnerClass import *
import re


def configureAndRunPostprocessing():
    runner = Runner()
    pathToRuns = os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-3-layers-2-subbands")
    # Adding '' at the end to terminate path with /
    directories = [
        os.path.join(pathToRuns, dir, "")
        for dir in os.listdir(pathToRuns)
        if re.match("RUN_E_Fermi_80.0.*", dir)
    ]
    enable_sc = ("sc_gap_calculation", "enable_sc_gap_calc", False)
    enable_chern = ("chern_number_calculation", "enable_chern_number_calc", False)
    enable_dos = ("dos_calculation", "enable_dos_calc", True)
    for dir in directories:
        nmlDirectorySC = ("sc_gap_calculation", "enable_dos_calc", dir)
        nmlDirectoryChern = (
            "chern_number_calculation",
            "path_to_run_dir_chern_number",
            dir,
        )
        nmlDirectoryDos = ("dos_calculation", "path_to_run_dir_dos", dir)
        runner.run_slurm_postprocessing(
            dir, [enable_dos, nmlDirectoryDos], True
        )


def configureAndRunSc():
    runner = Runner()

    # Setting number of sublattices and potential of layers
    n_sublattices = 3
    Sublat_param = ("discretization", "SUBLATTICES", n_sublattices)
    const_v_layer = 1.946e3
    V_layer_param = (
        "physical_params",
        "V_layer",
        [const_v_layer for _ in range(n_sublattices)],
    )
    n_subbands = 2
    Subband_param = ("discretization", "SUBBANDS", n_subbands)
    Subband_energies_param = ("physical_params", "Subband_energies", [0.0, 0.065e3])
    # Fermi energy
    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = 0.03e3
    Ef_max = 0.13e3
    Ef_steps = 10
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # J_SC
    nml_name = "physical_params"
    param_name = "J_SC"
    J_min = 0.15e3
    J_max = 0.29e3
    J_steps = 7
    dJ = abs(J_max - J_min) / J_steps
    J_table = [(nml_name, param_name, J_min + i * dJ) for i in range(J_steps + 1)]

    for Ef in Fermi_table:
        for J_sc in J_table:
            runner.run_slurm_param_value(
                paramValuePairs=[
                    Ef,
                    J_sc,
                    V_layer_param,
                    Sublat_param,
                    Subband_param,
                    Subband_energies_param,
                ],
                runsDir="KTO-SC/KTO-3-layers-2-subbands",
                material="KTO",
                isAres=True,
            )


def main():
    #configureAndRunSc()
    configureAndRunPostprocessing()

if __name__ == "__main__":
    main()
