import f90nml
import textwrap


class RunnerConfig:
    def __init__(self):

        # Using dedent to remove indentation and align every row with the first column of file
        self.job_header = textwrap.dedent(
            """\
        #!/bin/bash
        #SBATCH --job-name=LAO_STO              # Job name
        #SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes
        #SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node
        #SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)
        #SBATCH --mem-per-cpu=3800MB         # Memory (i.e. RAM) per processor
        #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory
        """
        )

        self.job_header_ares = textwrap.dedent(
            """\
        #!/bin/bash -l
        ## Job name
        #SBATCH -J KTO-SC
        ## Number of allocated nodes
        #SBATCH -N 1
        ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
        #SBATCH --ntasks-per-node=10
        ## Memory allocated per core (default is 5GB)
        #SBATCH --mem-per-cpu=3800MB
        ## Max task execution time (format is HH:MM:SS)
        #SBATCH --time=168:00:00
        ## Name of grant to which resource usage will be charged
        #SBATCH -A plglaosto111-cpu
        ## Name of partition
        #SBATCH -p plgrid-long
        ## Name of file to which standard output will be redirected
        #SBATCH --output="output.out"
        ## Name of file to which the standard error stream will be redirected
        #SBATCH --error="error.err"
        """
        )

    def LAO_STO_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&physical_params \
                T = 0.0, \
                t_D = 0.65e3, \
                t_I = 0.05e3, \
                t_Rashba = 0.004e3, \
                lambda_SOC = 0.265e3, \
                DELTA_TRI = -0.01e3, \
                v = 0.0e3, \
                V_pdp = 0.028e3, \
                V_pds = -0.065e3, \
                J_SC = 0.0e3, \
                J_SC_PRIME = 0.0e3, \
                J_SC_NNN = 0.0e3, \
                J_SC_PRIME_NNN = 0.0e3, \
                U_HUB = 0e3, \
                V_HUB = 0e3, \
                E_Fermi = -1.0e3 / \
            &discretization \
                k1_steps = 100, \
                k2_steps = 100 / \
            &self_consistency \
                read_gamma_from_file = .FALSE., \
                path_to_gamma_start = , \
                read_charge_from_file = .FALSE., \
                path_to_charge_start = , \
                gamma_start = 1., \
                gamma_nnn_start = 0., \
                charge_start = 0.1, \
                max_sc_iter = 100, \
                sc_alpha = 0.2, \
                sc_alpha_adapt = 1., \
                gamma_eps_convergence = 1e-4, \
                charge_eps_convergence = 1e-4 / \
            &romberg_integration \
                romb_eps_x = 1e-4, \
                interpolation_deg_x = 3, \
                max_grid_refinements_x = 14, \
                romb_eps_y = 1e-4, \
                interpolation_deg_y = 3, \
                max_grid_refinements_y = 14 /"
        )
        return params_nml

    def LAO_STO_postprocessing_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&sc_gap_calculation \
                enable_sc_gap_calc = .FALSE., \
                path_to_run_dir_sc_gap = , \
                dE_sc_gap = 1e-3, \
                Nk_points_sc_gap = 10000 / \
            &chern_number_calculation \
                enable_chern_number_calc = .FALSE., \
                path_to_run_dir_chern_number = , \
                Nk_points_chern_number = 15000 /\
            &dispersion_relation_calculation\
                enable_dispersion_relation_calc = .FALSE.,\
                path_to_run_dir_dispersion_relation = ,\
                Nk_points_dispersion_relation = 500 /\
            &dos_calculation\
                enable_dos_calc = .FALSE.,\
                path_to_run_dir_dos = ,\
                E_DOS_min = -1.7e3,\
                E_DOS_max = 1.7e3,\
                dE0 = 1.,\
                zeta_DOS = 1e-4,\
                Nk_points_dos = 2000 /"
        )
        return params_nml
