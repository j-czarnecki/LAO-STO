from RunnerClass import *
import re


def runTemperatureDependence():
    runner = Runner()
    pathToT0 = os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-fit-E_Fermi_100-150meV")

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
    Ef_min = 0.05e3
    Ef_max = 0.1e3
    Ef_steps = 5
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # J_SC
    nml_name = "physical_params"
    param_name = "J_SC"
    J_min = 0.36e3
    J_max = 0.37e3
    J_steps = 1
    dJ = abs(J_max - J_min) / J_steps
    J_table = [(nml_name, param_name, J_min + i * dJ) for i in range(J_steps + 1)]

    nml_name = "physical_params"
    param_name = "T"
    T_min = 0.25
    T_max = 6.0
    T_steps = 23
    dT = abs(T_max - T_min) / T_steps
    T_table = [(nml_name, param_name, T_min + i * dT) for i in range(T_steps + 1)]

    read_gamma = ("self_consistency", "read_gamma_from_file", True)
    read_charge = ("self_consistency", "read_charge_from_file", True)

    for T in T_table:
        for Ef in Fermi_table:
            for J_sc in J_table:
                t0DirGamma =os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_{J_sc[2]}_SUBLATTICES_3_SUBBANDS_2", "OutputData", "Gamma_SC_final.dat")
                t0DirCharge =os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_{J_sc[2]}_SUBLATTICES_3_SUBBANDS_2", "OutputData", "Charge_dens_final.dat")
                gammaDirSetting = ("self_consistency", "path_to_gamma_start", t0DirGamma)
                chargeDirSetting = ("self_consistency", "path_to_charge_start", t0DirCharge)

                runner.run_slurm_param_value(
                    paramValuePairs=[
                        T,
                        Ef,
                        J_sc,
                        V_layer_param,
                        Sublat_param,
                        Subband_param,
                        Subband_energies_param,
                        read_gamma,
                        read_charge,
                        gammaDirSetting,
                        chargeDirSetting
                    ],
                    runsDir="KTO-SC/KTO-fit-temperature",
                    material="KTO",
                    isAres=True,
                )

def configureAndRunPostprocessing():
    runner = Runner()
    pathToRuns = os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC")
    # Adding '' at the end to terminate path with /
    # directories = [
    #     os.path.join(pathToRuns, dir, "")
    #     for dir in os.listdir(pathToRuns)
    #     if re.match("RUN.*", dir)
    # ]
    directories = [
        os.path.join(pathToRuns, f"RUN_E_Fermi_{ef}.0_J_SC_110.0", "")
        for ef in [-1026, -966]]
    enable_sc = ("sc_gap_calculation", "enable_sc_gap_calc", True)
    enable_chern = ("chern_number_calculation", "enable_chern_number_calc", False)
    enable_dos = ("dos_calculation", "enable_dos_calc", False)
    for dir in directories:
        nmlDirectorySC = ("sc_gap_calculation", "path_to_run_dir_sc_gap", dir)
        nmlDirectoryChern = (
            "chern_number_calculation",
            "path_to_run_dir_chern_number",
            dir,
        )
        nmlDirectoryDos = ("dos_calculation", "path_to_run_dir_dos", dir)
        runner.run_slurm_postprocessing(
            dir, [enable_sc, nmlDirectorySC], True
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
    Ef_min = 0.05e3
    Ef_max = 0.15e3
    Ef_steps = 10
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # J_SC
    nml_name = "physical_params"
    param_name = "J_SC"
    J_min = 0.38e3
    J_max = 0.4e3
    J_steps = 2
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
                runsDir="KTO-SC/KTO-fit-E_Fermi_100-150meV",
                material="KTO",
                isAres=True,
            )


def main():
    #runTemperatureDependence()
    #configureAndRunSc()
    configureAndRunPostprocessing()

if __name__ == "__main__":
    main()
