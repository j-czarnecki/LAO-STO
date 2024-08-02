import f90nml


class RunnerConfig:
    def __init__(self):

        self.job_header = f"#!/bin/bash\n\
        #SBATCH --job-name=LAO_STO              # Job name\n\
        #SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes\n\
        #SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node\n\
        #SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)\n\
        #SBATCH --mem-per-cpu=3800MB         # Memory (i.e. RAM) per processor\n\
        #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory\n"



        self.job_header_ares = f'#!/bin/bash -l\n\
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

    def LAO_STO_default_nml(self):
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

    def LAO_STO_postprocessing_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(f'&sc_gap_calculation \
                                        enable_sc_gap_calc = .FALSE., \
                                        path_to_run_dir_sc_gap = , \
                                        dE_sc_gap = 1e-3, \
                                        Nk_points_sc_gap = 30000 / \
                                    &chern_number_calculation \
                                        enable_chern_number_calc = .FALSE., \
                                        path_to_run_dir_chern_number = , \
                                        Nk_points_chern_number = 30000 /')
        return params_nml