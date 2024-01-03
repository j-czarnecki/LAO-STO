import multiprocessing
import subprocess
import f90nml
import os
import sys


def LAO_STO_default_nml():
    parser = f90nml.Parser()
    params_nml = parser.reads(f'&physical_params \
                                    T = 0.00001, \
                                    t_D = 0.5e3, \
                                    t_I = 0.04e3, \
                                    lambda_SOC = 0.01e3, \
                                    DELTA_TRI = -0.005e3, \
                                    v = 0.2e3, \
                                    V_pdp = 0.028e3, \
                                    V_pds = -0.065e3, \
                                    J_SC = 0.175e3, \
                                    J_SC_PRIME = 0.0175e3, \
                                    U_HUB = 2e4, \
                                    V_HUB = 2e4, \
                                    E_Fermi = -0.5e3 / \
                                &discretization \
                                    k1_steps = 6000, \
                                    k2_steps = 6000 / \
                                &self_consistency \
                                    max_sc_iter = 1000, \
                                    sc_alpha = 0.2, \
                                    eps_convergence = 1e-7 /')
    return params_nml

def run_single_Ef(E_Fermi):
    print("Running for E_F = ", E_Fermi)
    path = f'RUN_Ef_' + str(E_Fermi)
    output_dir = f'OutputData'
    if not os.path.exists(path):
        os.mkdir(path)
        os.mkdir(os.path.join(path,output_dir))

    os.chdir(path)
    nml = LAO_STO_default_nml()
    nml['physical_params']['E_Fermi'] = E_Fermi
    nml['self_consistency']['max_sc_iter'] = 10
    with open('input.nml', 'w') as nml_file:
        f90nml.write(nml, nml_file, sort = False)
    simulate = subprocess.run(f'../LAO_STO.x')
    os.chdir(f'../')


def run_slurm(E_Fermi):
    job_header = f"#!/bin/bash\n\
#SBATCH --job-name=LAO_STO              # Job name\n\
#SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes\n\
#SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node\n\
#SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)\n\
#SBATCH --mem-per-cpu=8GB          # Memory (i.e. RAM) per processor\n\
#SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory\n"


    path = os.path.join(os.getcwd(), f'RUN_Ef_' + str(E_Fermi))
    output_dir = f'OutputData'
    if not os.path.exists(path):
        os.mkdir(path)
        os.mkdir(os.path.join(path,output_dir))
    os.chdir(path)
    nml = LAO_STO_default_nml()
    nml['physical_params']['E_Fermi'] = E_Fermi
    with open('input.nml', 'w') as nml_file:
        f90nml.write(nml, nml_file, sort = False)
    with open('job.sh', 'w') as job_file:
        print(job_header, file = job_file)
        print('cd ' + path, file = job_file)
        print(f'../LAO_STO.x', file = job_file)

    simulate = subprocess.run(["sbatch", "job.sh"])
    os.chdir(f'../')

if __name__ == '__main__':
    Ef_min = -1.5e3
    Ef_max = -0.5e3
    Ef_steps = 100
    dE = abs(Ef_max - Ef_min)/Ef_steps
    Fermi_table = [Ef_min + i*dE for i in range(Ef_steps + 1)]
    
    for Ef in Fermi_table:
        run_slurm(Ef)

    # try:
    #     with multiprocessing.Pool(processes=8) as pool:
    #         pool.map(run_single_Ef, Fermi_table)
    # except KeyboardInterrupt:
    #     print("Ctrl+C detected. Terminating processes.")
    #     pool.terminate()
    #     pool.join()
    #     print("Processes terminated.")
