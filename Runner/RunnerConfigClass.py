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
        #SBATCH -J v_tuning_NNN
        ## Number of allocated nodes
        #SBATCH -N 1
        ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
        #SBATCH --ntasks-per-node=48
        ## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
        ##SBATCH --mem-per-cpu=3800MB
        ## Memory allocated for whole job, comment if mem-per-cpu should be taken
        #SBATCH --mem=16GB
        ## Max task execution time (format is HH:MM:SS)
        #SBATCH --time=72:00:00
        ## Name of grant to which resource usage will be charged
        #SBATCH -A plgktosto111-cpu
        ## Name of partition
        #SBATCH -p plgrid
        ## Name of file to which standard output will be redirected
        #SBATCH --output="output.out"
        ## Name of file to which the standard error stream will be redirected
        #SBATCH --error="error.err"
        #module load GCC/13.2.0 OpenMPI/5.0.3 FlexiBLAS/3.3.1 ScaLAPACK/2.2.0-fb gimkl/2023b
        """
        )

    def LAO_STO_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&discretization \
                k1_steps = 100, \
                k2_steps = 100, \
                SUBLATTICES = 2, \
                SUBBANDS = 1 / \
            &physical_params \
                T = 0.0, \
                t_D = 0.5e3, \
                t_I = 0.04e3, \
                t_Rashba = 0.000e3, \
                lambda_SOC = 0.01e3, \
                delta_trigonal = -0.005e3, \
                zeta_tetragonal = 0.0, \
                orb_affected_tetragonal = 1, \
                v = 0.2e3, \
                V_pdp = 0.028e3, \
                V_pds = -0.065e3, \
                J_SC = 0.0e3, \
                J_SC_PRIME = 0.0e3, \
                J_SC_NNN = 0.0e3, \
                J_SC_PRIME_NNN = 0.0e3, \
                U_HUB = 0e3, \
                V_HUB = 0e3, \
                E_Fermi = -1.0e3, \
                V_layer = 1053.0, 1053.0, \
                Subband_energies = 0.0e3, \
                b_field = 0.0, 0.0, 0.0 / \
            &self_consistency \
                read_gamma_from_file = .FALSE., \
                path_to_gamma_start = '', \
                read_charge_from_file = .FALSE., \
                path_to_charge_start = '', \
                gamma_start = 1., \
                gamma_nnn_start = 1., \
                charge_start = 0.1, \
                max_sc_iter = 50, \
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

    def LAO_KTO_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&discretization \
                k1_steps = 100, \
                k2_steps = 100, \
                SUBLATTICES = 3, \
                SUBBANDS = 1 / \
            &physical_params \
                T = 0.0, \
                t_D = 0.65e3, \
                t_I = 0.05e3, \
                t_Rashba = 0.004e3, \
                lambda_SOC = 0.265e3, \
                delta_trigonal = 0.0e3, \
                zeta_tetragonal = 0.0, \
                orb_affected_tetragonal = 1, \
                v = 0.0e3, \
                V_pdp = 0.028e3, \
                V_pds = -0.065e3, \
                J_SC = 0.0e3, \
                J_SC_PRIME = 0.0e3, \
                J_SC_NNN = 0.0e3, \
                J_SC_PRIME_NNN = 0.0e3, \
                U_HUB = 0e3, \
                V_HUB = 0e3, \
                E_Fermi = -1.0e3, \
                V_layer = 1.946e3, 1.946e3, 1.946e3, \
                Subband_energies = 0.0e3, \
                b_field = 0.0, 0.0, 0.0 / \
            &self_consistency \
                read_gamma_from_file = .FALSE., \
                path_to_gamma_start = '', \
                read_charge_from_file = .FALSE., \
                path_to_charge_start = '', \
                gamma_start = 1., \
                gamma_nnn_start = 1., \
                charge_start = 0.1, \
                max_sc_iter = 30, \
                sc_alpha = 0.2, \
                sc_alpha_adapt = 1., \
                gamma_eps_convergence = 5e-3, \
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
                path_to_run_dir_sc_gap = '', \
                dE_sc_gap = 1e-3, \
                Nk_points_sc_gap = 4000, \
                Nk_points_sc_gap_refined = 4 / \
            &chern_number_calculation \
                enable_chern_number_calc = .FALSE., \
                path_to_run_dir_chern_number = '', \
                Nk_points_chern_number = 15000 /\
            &dispersion_relation_calculation\
                enable_dispersion_relation_calc = .FALSE.,\
                path_to_run_dir_dispersion_relation = '',\
                include_sc_in_dispersion = .TRUE.,\
                Nk_points_dispersion_relation = 500 /\
            &dos_calculation\
                enable_dos_calc = .FALSE.,\
                path_to_run_dir_dos = '',\
                E_DOS_min = -0.35,\
                E_DOS_max = 0.35,\
                dE0 = 2e-3,\
                zeta_DOS = 1.5e-2,\
                include_sc_in_dos = .TRUE.,\
                Nk_points_dos = 1000,\
                Nk_points_dos_refined = 10 /\
            &gamma_k_calculation \
                enable_gamma_k_calc = .FALSE., \
                path_to_run_dir_gamma_k = '/home/pwojcik/LAO-STO/', \
                Nk_points_gamma_k = 700 /"
    )
        return params_nml

    def DOS_fitter_yaml(self):
        fitConfig = {
            "runsDir": "KTO-SC/NicolasDosFitting/J2_0mT",
            "dosExpPath": "/net/people/plgrid/plgjczarnecki/NicolasFit/J2_exp_0mT.dat",
            "eMax": 0.25,
            "bounds":[
                [-0.3, 0.3], #gamma1
                [-0.3, 0.3], #gamma2
                [-0.3, 0.3] #gamma3
            ]
        }
        return fitConfig
