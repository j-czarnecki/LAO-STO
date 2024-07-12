import multiprocessing
import subprocess
import os
import sys
if os.path.exists('/net/home/pwojcik/.local/lib/python2.7/site-packages'):
    sys.path.insert(0, '/net/home/pwojcik/.local/lib/python2.7/site-packages')
import f90nml
import time

job_header = f"#!/bin/bash\n\
#SBATCH --job-name=LAO_STO              # Job name\n\
#SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes\n\
#SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node\n\
#SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)\n\
#SBATCH --mem-per-cpu=8GB          # Memory (i.e. RAM) per processor\n\
#SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory\n"

job_header_ares = f'#!/bin/bash -l\n\
## Job name\n\
#SBATCH -J LAO-STO\n\
## Number of allocated nodes\n\
#SBATCH -N 1\n\
## Number of tasks per node (by default this corresponds to the number of cores allocated per node)\n\
#SBATCH --ntasks-per-node=1\n\
## Memory allocated per core (default is 5GB)\n\
#SBATCH --mem-per-cpu=3800MB\n\
## Max task execution time (format is HH:MM:SS)\n\
#SBATCH --time=168:00:00\n\
## Name of grant to which resource usage will be charged\n\
#SBATCH -A plglaosto111-cpu\n\
## Name of partition\n\
#SBATCH -p plgrid-long\n\
## Name of file to which standard output will be redirected\n\
#SBATCH --output="output.out"\n\
## Name of file to which the standard error stream will be redirected\n\
#SBATCH --error="error.err"\n'


def LAO_STO_default_nml():
    parser = f90nml.Parser()
    params_nml = parser.reads(f'&physical_params \
                                    T = 0.0, \
                                    t_D = 0.5e3, \
                                    t_I = 0.04e3, \
                                    lambda_SOC = 0.01e3, \
                                    DELTA_TRI = -0.005e3, \
                                    v = 0.2e3, \
                                    V_pdp = 0.028e3, \
                                    V_pds = -0.065e3, \
                                    J_SC = 0.175e3, \
                                    J_SC_PRIME = 0.0175e3, \
                                    J_SC_NNN = 0.0e3, \
                                    J_SC_PRIME_NNN = 0.0e3, \
                                    U_HUB = 2e3, \
                                    V_HUB = 2e3, \
                                    E_Fermi = -0.1e3 / \
                                &discretization \
                                    k1_steps = 10, \
                                    k2_steps = 10 / \
                                &self_consistency \
                                    read_gamma_from_file = .FALSE., \
                                    path_to_gamma_start = , \
                                    read_charge_from_file = .FALSE., \
                                    path_to_charge_start = , \
                                    gamma_start = 1., \
                                    gamma_nnn_start = 0., \
                                    charge_start = 0.1, \
                                    max_sc_iter = 1000, \
                                    sc_alpha = 0.4, \
                                    sc_alpha_adapt = 1., \
                                    gamma_eps_convergence = 1e-4, \
                                    charge_eps_convergence = 1e-4 / \
                                &romberg_integration \
                                    romb_eps_x = 1e-6, \
                                    interpolation_deg_x = 4, \
                                    max_grid_refinements_x = 10, \
                                    romb_eps_y = 1e-6, \
                                    interpolation_deg_y = 4, \
                                    max_grid_refinements_y = 10 /')
    return params_nml

def run_slurm_param_value(paramValuePairs, isAres: bool = False):
    '''
    Sets all parameters given in key-value pairs.
    No need to specify J_SC_PRIME, since it is always set
    to be J_SC / 10
    '''
    if isAres:
        pathToAppend = f'/net/ascratch/people/plgjczarnecki/LAO-STO/RUN'
    else:
        pathToAppend = f'RUN'

    for pair in paramValuePairs:
        if pair[0] != 'self_consistency':
            pathToAppend = pathToAppend + f'_{pair[1]}_{pair[2]}'
    
    runner_cwd = os.getcwd()
    if isAres:
        path = pathToAppend
    else:
        path = os.path.join(runner_cwd, pathToAppend)

    output_dir = f'OutputData'
    if not os.path.exists(path):
        os.mkdir(path)
        os.mkdir(os.path.join(path,output_dir))
    os.chdir(path)

    nml = LAO_STO_default_nml() #creating default namelist
    for pair in paramValuePairs:
        nml[pair[0]][pair[1]] = pair[2] #editing all key-value pairs
    nml['physical_params']['J_SC_PRIME'] = nml['physical_params']['J_SC'] / 10. #No need to specify J_SC_PRIME
    nml['physical_params']['J_SC_PRIME_NNN'] = nml['physical_params']['J_SC_NNN'] / 10. #No need to specify J_SC_PRIME_NNN
    
    with open('input.nml', 'w') as nml_file:
        f90nml.write(nml, nml_file, sort = False)
    #setting up slurm script
    with open('job.sh', 'w') as job_file:
        if isAres:
            print(job_header_ares, file = job_file)
        else:
            print(job_header, file = job_file)

        print('cd ' + path, file = job_file)
        print(os.path.join(runner_cwd, 'LAO_STO.x'), file = job_file)

    #queue slurm job
    #simulate = subprocess.run(["sbatch", "job.sh"])
    os.chdir(runner_cwd)
    return path #for sequential runner

def run_sequential():
    '''
    Runs jobs sequentialy i.e. output of first job is the starting point for the second etc.
    '''
    nml_name = 'physical_params'
    param_name = 'E_Fermi'
    Ef_min = -1.1e3
    Ef_max = 0e3
    Ef_steps = 150
    dE = abs(Ef_max - Ef_min)/Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i*dE) for i in range(Ef_steps + 1)]

    isFirstIter = True
    previousPathRun = ''

    for Ef in Fermi_table:
        
        if isFirstIter:
            pathRun = run_slurm_param_value(
                [Ef,
                ('self_consistency', 'read_gamma_from_file', False),
                ('self_consistency', 'read_charge_from_file', False)],
                isAres=False)
            isFirstIter = False
        else:
            pathRun = run_slurm_param_value(
                [Ef,
                ('self_consistency', 'read_gamma_from_file', True),
                ('self_consistency', 'path_to_gamma_start', os.path.join(previousPathRun, 'OutputData/Gamma_SC_final.dat')),
                ('self_consistency', 'read_charge_from_file', True),
                ('self_consistency', 'path_to_charge_start', os.path.join(previousPathRun, 'OutputData/Charge_dens_final.dat'))],
                isAres=False)

        previousPathRun = pathRun

        while True:
            if os.path.exists(os.path.join(previousPathRun, 'OutputData/Gamma_SC_final.dat')):
                break
            else:
                time.sleep(30)

    print('!!! ALL SEQUENTIAL JOBS FINISHED !!!')

if __name__ == '__main__':

    run_sequential()

    #Fermi energy setting
    # nml_name = 'physical_params'
    # param_name = 'E_Fermi'
    # Ef_min = -1.05e3
    # Ef_max = -0.8e3
    # Ef_steps = 300
    # dE = abs(Ef_max - Ef_min)/Ef_steps
    # Fermi_table = [(nml_name, param_name, Ef_min + i*dE) for i in range(Ef_steps + 1)]
    # for Ef in Fermi_table:
    #     run_slurm_param_value([Ef], True)

    # try:
    #     with multiprocessing.Pool(processes=8) as pool:
    #         pool.map(run_single_Ef, Fermi_table)
    # except KeyboardInterrupt:
    #     print("Ctrl+C detected. Terminating processes.")
    #     pool.terminate()
    #     pool.join()
    #     print("Processes terminated.")

#def run_slurm(E_Fermi):
#     path = os.path.join(os.getcwd(), f'RUN_Ef_' + str(E_Fermi))
#     output_dir = f'OutputData'
#     if not os.path.exists(path):
#         os.mkdir(path)
#         os.mkdir(os.path.join(path,output_dir))
#     os.chdir(path)
#     nml = LAO_STO_default_nml()
#     nml['physical_params']['E_Fermi'] = E_Fermi
#     with open('input.nml', 'w') as nml_file:
#         f90nml.write(nml, nml_file, sort = False)
#     with open('job.sh', 'w') as job_file:
#         print(job_header, file = job_file)
#         print('cd ' + path, file = job_file)
#         print(f'../LAO_STO.x', file = job_file)

#     simulate = subprocess.run(["sbatch", "job.sh"])
#     os.chdir(f'../')  
          
#def run_single_Ef(E_Fermi):
#     print("Running for E_F = ", E_Fermi)
#     path = f'RUN_Ef_' + str(E_Fermi)
#     output_dir = f'OutputData'
#     if not os.path.exists(path):
#         os.mkdir(path)
#         os.mkdir(os.path.join(path,output_dir))

#     os.chdir(path)
#     nml = LAO_STO_default_nml()
#     nml['physical_params']['E_Fermi'] = E_Fermi
#     nml['self_consistency']['max_sc_iter'] = 10
#     with open('input.nml', 'w') as nml_file:
#         f90nml.write(nml, nml_file, sort = False)
#     simulate = subprocess.run(f'../LAO_STO.x')
#     os.chdir(f'../')