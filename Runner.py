import multiprocessing
import subprocess
import f90nml
import os
import sys


def LAO_STO_default_nml():
    parser = f90nml.Parser()
    params_nml = parser.reads(f'&physical_params \
                                    T = 0.01, \
                                    t_D = 0.5e3, \
                                    t_I = 0.04e3, \
                                    lambda_SOC = 0.01e3, \
                                    DELTA_TRI = -0.005e3, \
                                    v = 0.2e3, \
                                    V_pdp = 0.028e3, \
                                    V_pds = -0.065e3, \
                                    J_SC = 0.165e3, \
                                    J_SC_PRIME = 0.0165e3, \
                                    U_HUB = 2e4, \
                                    V_HUB = 2e4, \
                                    E_Fermi = -0.5e3 / \
                                &discretization \
                                    k1_steps = 5000, \
                                    k2_steps = 5000 / \
                                &self_consistency \
                                    max_sc_iter = 1000, \
                                    sc_alpha = 0.1, \
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


if __name__ == '__main__':
    Ef_min = -0.7e3
    Ef_max = -0.4e3
    Ef_steps = 30
    dE = abs(Ef_max - Ef_min)/Ef_steps
    Fermi_table = [Ef_min + i*dE for i in range(Ef_steps + 1)]
    
    with multiprocessing.Pool(processes=8) as pool:
        pool.map(run_single_Ef, Fermi_table)