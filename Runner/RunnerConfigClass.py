# This file is part of LAO-STO.
#
# Copyright (C) 2025 Julian Czarnecki
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# If you use this code for scientific research, please cite:
# J. Czarnecki et. al.,
# "Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface",
# arXiv:2508.05075 (2025).
# https://arxiv.org/abs/2508.05075

import f90nml
import textwrap
from collections import OrderedDict

class RunnerConfig:
    def __init__(self):

        # Using dedent to remove indentation and align every row with the first column of file
        self.job_header = {
            "default": textwrap.dedent(
                """\
                #!/bin/bash
                ##### Amount of cores per task
                #SBATCH --cpus-per-task=1
                ##### Partition name
                #SBATCH -p cpu
                ##### Name of job in queuing system
                #SBATCH --job-name=KTO-B_planar
                #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory
                """
                ),
            "ares": textwrap.dedent(
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
                """
                ),
            "helios": textwrap.dedent(
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
                module load GCC/13.2.0 OpenMPI/5.0.3 FlexiBLAS/3.3.1 ScaLAPACK/2.2.0-fb gimkl/2023b
                """
                ),
        }

        self.defaultSelfConsistencyNmlDict = OrderedDict([
            ("discretization", OrderedDict([
                ("k1_steps", 100),
                ("k2_steps", 100),
                ("SUBLATTICES", 2),
                ("SUBBANDS", 1),
              ])),
            ("physical_params", OrderedDict([
                ("T", 0.0),
                ("t_D", 0.5e3),
                ("t_I", 0.04e3),
                ("t_Rashba", 0.000e3),
                ("lambda_SOC", 0.01e3),
                ("delta_trigonal", -0.005e3),
                ("zeta_tetragonal", 0.0e3),
                ("orb_affected_tetragonal", 1),
                ("v", 0.0e3),
                ("V_pdp", 0.028e3),
                ("V_pds", -0.065e3),
                ("J_SC", 0.0e3),
                ("J_SC_PRIME", 0.0e3),
                ("J_SC_NNN", 0.0e3),
                ("J_SC_PRIME_NNN", 0.0e3),
                ("U_HUB", 0.0e3),
                ("V_HUB", 0.0e3),
                ("E_Fermi", 0.0e3),
                ("V_layer", [0.0, 0.0]),
                ("Subband_energies", 0.0e3),
                ("B_magnitude", 0.0),
                ("B_theta", 90.0),
                ("B_phi", 0.0),
            ])),
            ("self_consistency", OrderedDict([
                ("read_gamma_from_file", False),
                ("path_to_gamma_start", ''),
                ("read_charge_from_file", False),
                ("path_to_charge_start", ''),
                ("gamma_start", 1.),
                ("gamma_nnn_start", 1.),
                ("charge_start", 0.1),
                ("max_sc_iter", 50),
                ("sc_alpha", 0.2),
                ("sc_alpha_adapt", 1.),
                ("gamma_eps_convergence", 1e-4),
                ("charge_eps_convergence", 1e-4),
            ])),
            ("romberg_integration", OrderedDict([
                ("romb_eps_x", 5e-4),
                ("interpolation_deg_x", 3),
                ("max_grid_refinements_x", 11),
                ("romb_eps_y", 5e-4),
                ("interpolation_deg_y", 3),
                ("max_grid_refinements_y", 11),
            ])),
          ])

        self.defaultPostNmlDict = OrderedDict([
            ("sc_gap_calculation", OrderedDict([
                ("enable_sc_gap_calc", False),
                ("path_to_run_dir_sc_gap", ""),
                ("dE_sc_gap", 1e-3),
                ("Nk_points_sc_gap", 8000),
                ("Nk_points_sc_gap_refined", 20),
            ])),
            ("chern_number_calculation", OrderedDict([
                ("enable_chern_number_calc", False),
                ("path_to_run_dir_chern_number", ""),
                ("Nk_points_chern_number", 15000),
            ])),
            ("dispersion_relation_calculation", OrderedDict([
                ("enable_dispersion_relation_calc", False),
                ("path_to_run_dir_dispersion_relation", ""),
                ("include_sc_in_dispersion", True),
                ("Nk_points_dispersion_relation", 500),
            ])),
            ("dos_calculation", OrderedDict([
                ("enable_dos_calc", False),
                ("path_to_run_dir_dos", ""),
                ("E_DOS_min", -0.35),
                ("E_DOS_max", 0.35),
                ("dE0", 2e-3),
                ("zeta_DOS", 1.5e-2),
                ("include_sc_in_dos", True),
                ("Nk_points_dos", 1000),
                ("Nk_points_dos_refined", 10),
            ])),
            ("gamma_k_calculation", OrderedDict([
                ("enable_gamma_k_calc", False),
                ("path_to_run_dir_gamma_k", "/home/pwojcik/LAO-STO/"),
                ("Nk_points_gamma_k", 700),
            ])),
            ("projections_calculation", OrderedDict([
                ("enable_projections_calc", False),
                ("path_to_run_dir_projections", "/home/pwojcik/LAO-STO/"),
                ("Nr_points_projections", 50),
                ("Nphi_points_projections", 50),
            ])),
        ])

    def LAO_STO_default_nml(self):
        nml = f90nml.Namelist(self.defaultSelfConsistencyNmlDict)
        nml['physical_params']['t_D'] = 0.5e3
        nml['physical_params']['t_I'] = 0.04e3
        nml['physical_params']['lambda_SOC'] = 0.01e3
        nml['physical_params']['delta_trigonal'] = -0.005e3
        nml['physical_params']['v'] = 0.2e3
        nml['physical_params']['V_layer'] = [1053.0, 1053.0]
        nml['romberg_integration']['max_grid_refinements_x'] = 12
        nml['romberg_integration']['max_grid_refinements_y'] = 12
        return nml

    def LAO_KTO_default_nml(self):
        nml = f90nml.Namelist(self.defaultSelfConsistencyNmlDict)
        nml['discretization']['SUBLATTICES'] = 3
        nml['physical_params']['t_D'] = 0.65e3
        nml['physical_params']['t_I'] = 0.05e3
        nml['physical_params']['t_Rashba'] = 0.004e3
        nml['physical_params']['lambda_SOC'] = 0.265e3
        nml['physical_params']['delta_trigonal'] = 0.01e3
        nml['physical_params']['V_layer'] = [1946.0, 1946.0, 1946.0]
        return nml

    def LAO_STO_postprocessing_default_nml(self):
        return f90nml.Namelist(self.defaultPostNmlDict)

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
