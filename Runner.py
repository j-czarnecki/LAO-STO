import multiprocessing
import subprocess
import os
import sys
if os.path.exists('/net/home/pwojcik/.local/lib/python2.7/site-packages'):
    sys.path.insert(0, '/net/home/pwojcik/.local/lib/python2.7/site-packages')
import f90nml

job_header = f"#!/bin/bash\n\
#SBATCH --job-name=LAO_STO              # Job name\n\
#SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes\n\
#SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node\n\
#SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)\n\
#SBATCH --mem-per-cpu=8GB          # Memory (i.e. RAM) per processor\n\
#SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory\n"


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
                                    U_HUB = 2e4, \
                                    V_HUB = 2e4, \
                                    E_Fermi = -0.1e3 / \
                                &discretization \
                                    k1_steps = 10, \
                                    k2_steps = 10 / \
                                &self_consistency \
                                    gamma_start = 1., \
                                    max_sc_iter = 1000, \
                                    sc_alpha = 0.4, \
                                    eps_convergence = 1e-4 / \
                                &romberg_integration \
                                    romb_eps_x = 1e-6, \
                                    interpolation_deg_x = 4, \
                                    max_grid_refinements_x = 10, \
                                    romb_eps_y = 1e-6, \
                                    interpolation_deg_y = 4, \
                                    max_grid_refinements_y = 10 /')
    return params_nml

def run_slurm_param_value(paramValuePairs):
    '''
    Sets all parameters given in key-value pairs.
    No need to specify J_SC_PRIME, since it is always set
    to be J_SC / 10
    '''

    pathToAppend = f'RUN'
    for pair in paramValuePairs:
        pathToAppend = pathToAppend + f'_{pair[1]}_{pair[2]}'
    path = os.path.join(os.getcwd(), pathToAppend)

    output_dir = f'OutputData'
    if not os.path.exists(path):
        os.mkdir(path)
        os.mkdir(os.path.join(path,output_dir))
    os.chdir(path)

    nml = LAO_STO_default_nml() #creating default namelist
    for pair in paramValuePairs:
        nml[pair[0]][pair[1]] = pair[2] #editing all key-value pairs
    nml['physical_params']['J_SC_PRIME'] = nml['physical_params']['J_SC'] / 10. #No need to specify J_SC_PRIME

    
    with open('input.nml', 'w') as nml_file:
        f90nml.write(nml, nml_file, sort = False)
    #setting up slurm script
    with open('job.sh', 'w') as job_file:
        print(job_header, file = job_file)
        print('cd ' + path, file = job_file)
        print(f'../LAO_STO.x', file = job_file)

    #queue slurm job
    simulate = subprocess.run(["sbatch", "job.sh"])
    os.chdir(f'../')

if __name__ == '__main__':

    #Fermi energy setting
    nml_name = 'physical_params'
    param_name = 'E_Fermi'
    Ef_min = -1.1e3
    Ef_max = -1.0e3
    Ef_steps = 4
    dE = abs(Ef_max - Ef_min)/Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i*dE) for i in range(Ef_steps + 1)]

    #J_SC setting
    nml_name = 'physical_params'
    param_name = 'J_SC'
    J_SC_min = 0.
    J_SC_max = 0.2e3
    J_SC_steps = 4
    dJ = abs(J_SC_max - J_SC_min)/J_SC_steps
    J_SC_table = [(nml_name, param_name, J_SC_min + i*dJ) for i in range(J_SC_steps + 1)]
    
    for Ef in Fermi_table:
        for J_SC in J_SC_table:
            paramsList = [Ef] + [J_SC]
            run_slurm_param_value(paramsList)

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