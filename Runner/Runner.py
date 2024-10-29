from RunnerClass import *
import re


def main():
    runner = Runner()

    pathToRuns = "/net/ascratch/people/plgjczarnecki/LAO-STO-NNN-test/"
    #Adding '' at the end to terminate path with /
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

    # runner.run_sequential()

    # #Run as a function of temperature
    # # Fermi energy setting
    # nml_name = "physical_params"
    # param_name = "T"
    # T_min = 0.73
    # T_max = 1.49
    # T_steps = 38
    # dT = abs(T_max - T_min) / T_steps
    # T_table = [(nml_name, param_name, T_min + i * dT) for i in range(T_steps + 1)]  

    # # # Fermi energy setting
    # nml_name = "physical_params"
    # param_name = "E_Fermi"
    # Ef_min = -0.972e3
    # Ef_max = -0.964e3
    # Ef_steps = 2
    # dE = abs(Ef_max - Ef_min) / Ef_steps
    # Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # nml_name = "physical_params"
    # param_name = 'J_SC'
    # J_SC = 0.11e3
    # J_SC_nml = (nml_name, param_name, J_SC)

    # T0_path_prefix = '/net/ascratch/people/plgjczarnecki/LAO-STO-E_Fermi_J_SC'
    # for Ef in Fermi_table:
    #     for T in T_table:
    #         T0_path = os.path.join(T0_path_prefix, f'RUN_E_Fermi_{Ef[2]}_J_SC_110.0')
    #         gamma_path = os.path.join(T0_path, 'OutputData', 'Gamma_SC_iter.dat')
    #         charge_path = os.path.join(T0_path, 'OutputData', 'Charge_dens_iter.dat')
    #         gamma_to_nml = ('self_consistency', 'path_to_gamma_start', gamma_path)
    #         charge_to_nml = ('self_consistency', 'path_to_charge_start', charge_path)
    #         gamma_from_path_true = ('self_consistency', 'read_gamma_from_file', True)
    #         charge_from_path_true = ('self_consistency', 'read_charge_from_file', True)
    #         runner.run_slurm_param_value([Ef,
    #                                       J_SC_nml,
    #                                       T,
    #                                       gamma_from_path_true,
    #                                       gamma_to_nml,
    #                                       charge_from_path_true,
    #                                       charge_to_nml],
    #                                       isAres=True)


    # # for Ef in Fermi_table:
    # #     runner.run_slurm_param_value([Ef], True)
    #Fermi energy
    # nml_name = "physical_params"
    # param_name = "E_Fermi"
    # Ef_min = -0.972e3
    # Ef_max = -0.964e3
    # Ef_steps = 2
    # dE = abs(Ef_max - Ef_min) / Ef_steps
    # Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]


    # nml_name = "physical_params"
    # param_name = "J_SC_NNN"
    # J_min = 0.10e3
    # J_max = 0.15e3
    # J_steps = 1
    # dJ = abs(J_max - J_min) / J_steps
    # J_table = [(nml_name, param_name, J_min + i * dJ) for i in range(J_steps + 1)]
    # for Ef in Fermi_table:
    #     for J_sc in J_table:
    #         runner.run_slurm_param_value([Ef, J_sc], True)

    # try:
    #     with multiprocessing.Pool(processes=8) as pool:
    #         pool.map(run_single_Ef, Fermi_table)
    # except KeyboardInterrupt:
    #     print("Ctrl+C detected. Terminating processes.")
    #     pool.terminate()
    #     pool.join()
    #     print("Processes terminated.")


if __name__ == "__main__":
    main()
