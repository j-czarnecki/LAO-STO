from RunnerClass import *
import re

def main():
    runner = Runner()

    
    pathToRuns = '/net/ascratch/people/plgjczarnecki/LAO-STO-v0'
    #Adding '' at the end to terminate path with /
    directories = [os.path.join(pathToRuns, dir, '') for dir in os.listdir(pathToRuns) if re.match('RUN_.*', dir)]
    enable_sc = ('sc_gap_calculation', 'enable_sc_gap_calc', True)
    enable_chern = ('chern_number_calculation', 'enable_chern_number_calc', True)
    for dir in directories:
        nmlDirectorySC = ('sc_gap_calculation', 'path_to_run_dir_sc_gap', dir)
        nmlDirectoryChern = ('chern_number_calculation', 'path_to_run_dir_chern_number', dir)
        runner.run_slurm_postprocessing(dir, [enable_sc, nmlDirectorySC, enable_chern, nmlDirectoryChern])


    #runner.run_sequential()

    #Fermi energy setting
    # nml_name = 'physical_params'
    # param_name = 'E_Fermi'
    # Ef_min = -1.05e3
    # Ef_max = -0.8e3
    # Ef_steps = 300
    # dE = abs(Ef_max - Ef_min)/Ef_steps
    # Fermi_table = [(nml_name, param_name, Ef_min + i*dE) for i in range(Ef_steps + 1)]
    # for Ef in Fermi_table:
    #     runner.run_slurm_param_value([Ef], True)

    # try:
    #     with multiprocessing.Pool(processes=8) as pool:
    #         pool.map(run_single_Ef, Fermi_table)
    # except KeyboardInterrupt:
    #     print("Ctrl+C detected. Terminating processes.")
    #     pool.terminate()
    #     pool.join()
    #     print("Processes terminated.")

if __name__ == '__main__':
    main()