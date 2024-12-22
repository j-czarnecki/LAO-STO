from RunnerClass import *
import re


def configureAndRunPostprocessing():
    runner = Runner()
    pathToRuns = "/net/ascratch/people/plgjczarnecki/LAO-STO-test/"
    # Adding '' at the end to terminate path with /
    directories = [
        os.path.join(pathToRuns, dir, "")
        for dir in os.listdir(pathToRuns)
        if re.match("RUN_.*", dir)
    ]
    enable_sc = ("sc_gap_calculation", "enable_sc_gap_calc", True)
    enable_chern = ("chern_number_calculation", "enable_chern_number_calc", False)
    for dir in directories:
        nmlDirectorySC = ("sc_gap_calculation", "path_to_run_dir_sc_gap", dir)
        nmlDirectoryChern = (
            "chern_number_calculation",
            "path_to_run_dir_chern_number",
            dir,
        )
        runner.run_slurm_postprocessing(
            dir, [enable_sc, nmlDirectorySC, enable_chern, nmlDirectoryChern], True
        )


def configureAndRunSc():
    runner = Runner()
    # Fermi energy
    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = -1.405e3
    Ef_max = -1.25e3
    Ef_steps = 31
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # J_SC
    nml_name = "physical_params"
    param_name = "J_SC"
    J_min = 0.1e3
    J_max = 0.4e3
    J_steps = 3
    dJ = abs(J_max - J_min) / J_steps
    J_table = [(nml_name, param_name, J_min + i * dJ) for i in range(J_steps + 1)]
    for Ef in Fermi_table:
        for J_sc in J_table:
            runner.run_slurm_param_value(
                paramValuePairs=[Ef, J_sc],
                rusnDir="KTO-test",
                material="KTO",
                isAres=True,
            )


def main():
    configureAndRunSc()
    # configureAndRunPostprocessing()


if __name__ == "__main__":
    main()
