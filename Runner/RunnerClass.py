import multiprocessing
import subprocess
import os
import sys

if os.path.exists("/net/home/pwojcik/.local/lib/python2.7/site-packages"):
    sys.path.insert(0, "/net/home/pwojcik/.local/lib/python2.7/site-packages")
import time
from RunnerConfigClass import *
import yaml

SCRATCH_PATH = os.getenv('SCRATCH')
HOME_PATH = os.getenv('HOME')

class Runner(RunnerConfig):
    def __init__(self):
        RunnerConfig.__init__(self)

    def run_slurm_param_value(
        self, paramValuePairs: list, runsDir: str, material: str, machine: str = "default"
    ):
        """
        Sets all parameters given in key-value pairs.
        No need to specify J_SC_PRIME, since it is always set
        to be J_SC / 10.
        Also starting values of gamma parameters are automatically set to 0 if corresponding J is 0.
        """
        pathToAppend = os.path.join(SCRATCH_PATH, runsDir)
        os.makedirs(pathToAppend, exist_ok=True)
        pathToAppend = os.path.join(pathToAppend, "RUN")

        for pair in paramValuePairs:
            if pair[0] != "self_consistency" and not isinstance(pair[2], list):
                pathToAppend = pathToAppend + f"_{pair[1]}_{pair[2]}"

        runner_cwd = os.getcwd()
        path = pathToAppend

        output_dir = f"OutputData"
        if not os.path.exists(path):
            os.mkdir(path)
            os.mkdir(os.path.join(path, output_dir))
        os.chdir(path)

        # Getting namelist with parameters
        nml = f90nml.Namelist()
        if material == "STO":
            nml = self.LAO_STO_default_nml()
        elif material == "KTO":
            nml = self.LAO_KTO_default_nml()

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        self.__set_derived_values(nml)

        with open("input.nml", "w") as nml_file:
            f90nml.write(nml, nml_file, sort=False)
        # setting up slurm script
        with open("job.sh", "w") as job_file:
            print(self.job_header[machine], file=job_file)
            print("cd " + path, file=job_file)
            print(os.path.join(runner_cwd, "..", "bin", "LAO_STO.x"), file=job_file)

        # queue slurm job
        simulate = subprocess.run(["sbatch", "job.sh"])
        os.chdir(runner_cwd)
        return path  # for sequential runner

    def run_slurm_postprocessing(self, runDir, paramValuePairs, machine: str = "default"):

        runner_cwd = os.getcwd()
        os.chdir(runDir)

        nml = self.LAO_STO_postprocessing_default_nml()

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        with open("postprocessing_input.nml", "w") as nml_file:
            f90nml.write(nml, nml_file, sort=False)

        # setting up slurm script
        os.chdir(runDir)
        with open("job.sh", "w") as job_file:
            print(self.job_header[machine], file=job_file)

            print("cd " + runDir, file=job_file)
            print(
                os.path.join(runner_cwd, "..", "bin", "POST_LAO_STO.x"),
                file=job_file,
            )

        # queue slurm job
        simulate = subprocess.run(["sbatch", "job.sh"])
        os.chdir(runner_cwd)

    def run_slurm_dos_fitter(self,
                             fitConfig: dict,
                             paramValuePairs: list,
                             material: str,
                             machine: str = "default"):

        print(fitConfig)

        # Prepare directory structure
        os.makedirs(os.path.join(SCRATCH_PATH, fitConfig['runsDir']), exist_ok=True)
        os.makedirs(os.path.join(SCRATCH_PATH, fitConfig['runsDir'], "OutputData"), exist_ok=True)
        os.makedirs(os.path.join(SCRATCH_PATH, fitConfig['runsDir'], "Plots"), exist_ok=True)
        os.makedirs(os.path.join(SCRATCH_PATH, fitConfig['runsDir'], "DOS_train"), exist_ok=True)
        #A dummy charge dens file has to exists for postprocessing
        subprocess.run(f"cp {os.path.join(HOME_PATH,'LAO-STO', 'OutputData', 'Charge_dens_final.dat')} {os.path.join(SCRATCH_PATH, fitConfig['runsDir'], 'OutputData')}", shell=True, check=True)

        runnerCwd = os.getcwd()
        os.chdir(fitConfig["runsDir"])



        self.__create_and_write_input_nml(paramValuePairs, material)

        # Setup postprocessing_input
        nml = self.LAO_STO_postprocessing_default_nml()
        nml['dos_calculation']['enable_dos_calc'] = True
        nml['dos_calculation']['path_to_run_dir_dos'] = f"{fitConfig['runsDir']}/"
        with open("postprocessing_input.nml", "w") as nml_file:
            f90nml.write(nml, nml_file, sort=False)

        # Save fitConfig.yaml
        defaultFitConfig = self.DOS_fitter_yaml()
        for key in fitConfig:
            defaultFitConfig[key] = fitConfig[key]
        with open(os.path.join(SCRATCH_PATH, fitConfig["runsDir"], "fitConfig.yaml"), "w") as f:
            yaml.dump(defaultFitConfig, f, sort_keys=True)


        # setting up slurm script
        os.chdir(fitConfig["runsDir"])
        with open("job.sh", "w") as job_file:
            print(self.job_header[machine], file=job_file)

            print(f"cd {os.path.join(HOME_PATH, 'LAO-STO')}", file=job_file)
            print(f"python3 -m Fitter.mainFitter --config {os.path.join(fitConfig['runsDir'], 'fitConfig.yaml')}", file=job_file)

        # queue slurm job
        subprocess.run(["sbatch", "job.sh"], check=True)
        os.chdir(runnerCwd)

    def run_sequential(self):
        """
        ============= DEPRECATED =====================
        Runs jobs sequentialy i.e. output of first job is the starting point for the second etc.
        """
        nml_name = "physical_params"
        param_name = "E_Fermi"
        Ef_min = -1.1e3
        Ef_max = 0e3
        Ef_steps = 150
        dE = abs(Ef_max - Ef_min) / Ef_steps
        Fermi_table = [
            (nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)
        ]

        isFirstIter = True
        previousPathRun = ""

        for Ef in Fermi_table:

            # if isFirstIter:
            #     pathRun = self.run_slurm_param_value(
            #         [
            #             Ef,
            #             ("self_consistency", "read_gamma_from_file", False),
            #             ("self_consistency", "read_charge_from_file", False),
            #         ],
            #         isAres=False,
            #     )
            #     isFirstIter = False
            # else:
            #     pathRun = self.run_slurm_param_value(
            #         [
            #             Ef,
            #             ("self_consistency", "read_gamma_from_file", True),
            #             (
            #                 "self_consistency",
            #                 "path_to_gamma_start",
            #                 os.path.join(
            #                     previousPathRun, "OutputData/Gamma_SC_final.dat"
            #                 ),
            #             ),
            #             ("self_consistency", "read_charge_from_file", True),
            #             (
            #                 "self_consistency",
            #                 "path_to_charge_start",
            #                 os.path.join(
            #                     previousPathRun, "OutputData/Charge_dens_final.dat"
            #                 ),
            #             ),
            #         ],
            #         isAres=False,
            #     )

            previousPathRun = pathRun

            while True:
                if os.path.exists(
                    os.path.join(previousPathRun, "OutputData/Gamma_SC_final.dat")
                ):
                    break
                else:
                    time.sleep(30)

        print("!!! ALL SEQUENTIAL JOBS FINISHED !!!")

    def __set_derived_values(self, nml: f90nml.Namelist):
        """
        This method sets all values derived from other namelist parameters
        """
        # Setting inter-orbital pairing
        nml["physical_params"]["J_SC_PRIME"] = (
            nml["physical_params"]["J_SC"] / 10.0
        )  # No need to specify J_SC_PRIME
        nml["physical_params"]["J_SC_PRIME_NNN"] = (
            nml["physical_params"]["J_SC_NNN"] / 10.0
        )  # No need to specify J_SC_PRIME_NNN

        # Setting appropriate starting values
        if nml["physical_params"]["J_SC"] == 0.0:
            nml["self_consistency"]["gamma_start"] = 0.0
        if nml["physical_params"]["J_SC_NNN"] == 0.0:
            nml["self_consistency"]["gamma_nnn_start"] = 0.0

    def __create_and_write_input_nml(self, paramValuePairs: list, material: str) -> None:
        # Getting namelist with parameters
        nml = f90nml.Namelist()
        if material == "STO":
            nml = self.LAO_STO_default_nml()
        elif material == "KTO":
            nml = self.LAO_KTO_default_nml()

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        self.__set_derived_values(nml)

        with open("input.nml", "w") as nml_file:
            f90nml.write(nml, nml_file, sort=False)
