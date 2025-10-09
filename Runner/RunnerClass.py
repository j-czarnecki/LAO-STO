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

import multiprocessing
import subprocess
import os
import sys

if os.path.exists("/net/home/pwojcik/.local/lib/python2.7/site-packages"):
    sys.path.insert(0, "/net/home/pwojcik/.local/lib/python2.7/site-packages")
import time
from RunnerConfigClass import *
from OutputMocker import *
import yaml

SCRATCH_PATH = os.getenv('SCRATCH')
HOME_PATH = os.getenv('HOME')

class Runner(RunnerConfig):
    def __init__(self):
        RunnerConfig.__init__(self)

    def runSlurmParamValue(
        self, paramValuePairs: list[tuple[str, str, float|list[float]]], runsDir: str, material: str, machine: str = "default"
    ):
        """
        Sets all parameters given in key-value pairs.
        No need to specify J_SC_PRIME, since it is always set
        to be J_SC / 10.
        Also starting values of gamma parameters are automatically set to 0 if corresponding J is 0.
        """
        newRunPath = self.__createRunDirStructure(runsDir, paramValuePairs)
        runnerCwd = os.getcwd()
        os.chdir(newRunPath)

        # Getting namelist with parameters
        nml = self.__getMaterialNml(material)

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs


        with open("input.nml", "w") as nmlFile:
            f90nml.write(nml, nmlFile, sort=False)
        # setting up slurm script
        with open("job.sh", "w") as jobFile:
            print(self.jobHeader[machine], file=jobFile)
            print("cd " + newRunPath, file=jobFile)
            print(os.path.join(runnerCwd, "..", "bin", "LAO_STO.x"), file=jobFile)

        # queue slurm job
        simulate = subprocess.run(["sbatch", "job.sh"])
        os.chdir(runnerCwd)
        return newRunPath  # for sequential runner

    def runSlurmPostprocessing(self, runDir, paramValuePairs, machine: str = "default"):

        runnerCwd = os.getcwd()
        os.chdir(runDir)

        nml = self.getPostprocessingDefaultNml()

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        with open("postprocessing_input.nml", "w") as nmlFile:
            f90nml.write(nml, nmlFile, sort=False)

        # setting up slurm script
        os.chdir(runDir)
        with open("job.sh", "w") as jobFile:
            print(self.jobHeader[machine], file=jobFile)

            print("cd " + runDir, file=jobFile)
            print(
                os.path.join(runnerCwd, "..", "bin", "POST_LAO_STO.x"),
                file=jobFile,
            )

        # queue slurm job
        simulate = subprocess.run(["sbatch", "job.sh"])
        os.chdir(runnerCwd)

    def runSlurmDosFitter(self,
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



        self.__createAndWriteInputNml(paramValuePairs, material)

        # Setup postprocessing_input
        nml = self.getPostprocessingDefaultNml()
        nml['dos_calculation']['enable_dos_calc'] = True
        nml['dos_calculation']['path_to_run_dir_dos'] = f"{fitConfig['runsDir']}/"
        with open("postprocessing_input.nml", "w") as nmlFile:
            f90nml.write(nml, nmlFile, sort=False)

        # Save fitConfig.yaml
        defaultFitConfig = self.getDosFitterYamlDict()
        for key in fitConfig:
            defaultFitConfig[key] = fitConfig[key]
        with open(os.path.join(SCRATCH_PATH, fitConfig["runsDir"], "fitConfig.yaml"), "w") as f:
            yaml.dump(defaultFitConfig, f, sort_keys=True)


        # setting up slurm script
        os.chdir(fitConfig["runsDir"])
        with open("job.sh", "w") as jobFile:
            print(self.jobHeader[machine], file=jobFile)

            print(f"cd {os.path.join(HOME_PATH, 'LAO-STO')}", file=jobFile)
            print(f"python3 -m Fitter.mainFitter --config {os.path.join(fitConfig['runsDir'], 'fitConfig.yaml')}", file=jobFile)

        # queue slurm job
        subprocess.run(["sbatch", "job.sh"], check=True)
        os.chdir(runnerCwd)

    def runSlurmMockedOutputPostprocessing(self,
                                           paramValuePairs: list[tuple[str, str, float]],
                                           paramValuePairsPost: list[tuple[str, str, float]],
                                           gammaAmplitudesDict: dict[str, np.complex128],
                                           symmetriesWeightsDict: dict[str, dict[str, float]],
                                           runsDir: str,
                                           material: str,
                                           machine: str = "default"):
        """
        This method generates an artificial output - Gamma_SC_final.dat and Charge_dens_final.dat - in a given directory,
        along with desired input.nml and postrocessing_input.nml files. Eventually it runs postprocessing with given data.
        """
        # Creating new directory and input.nml
        newRunPath = self.__createRunDirStructure(runsDir, paramValuePairs)
        runnerCwd = os.getcwd()
        os.chdir(newRunPath)

        # Getting namelist with parameters
        nml = self.__getMaterialNml(material)

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        # Save input.nml
        with open("input.nml", "w") as nmlFile:
            f90nml.write(nml, nmlFile, sort=False)

        # Creating postprocessing_input.nml
        nmlPost = self.getPostprocessingDefaultNml()
        paramValuePairsPost.append(("sc_gap_calculation", "enable_sc_gap_calc", True))
        paramValuePairsPost.append(("sc_gap_calculation", "path_to_run_dir_sc_gap", f"{newRunPath}/"))
        for pair in paramValuePairsPost:
            nmlPost[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        with open("postprocessing_input.nml", "w") as nmlFile:
            f90nml.write(nmlPost, nmlFile, sort=False)

        # Setting up slurm script
        with open("job.sh", "w") as jobFile:
            print(self.jobHeader[machine], file=jobFile)

            print("cd " + newRunPath, file=jobFile)
            print(
                os.path.join(runnerCwd, "..", "bin", "POST_LAO_STO.x"),
                file=jobFile,
            )

        # Creating mocked output to be able to run postprocessing
        outputMocker = OutputMocker(newRunPath,
                                    nOrbs=3,
                                    nBands=nml['discretization']['subbands'],
                                    nSublats=nml['discretization']['sublattices'])
        outputMocker.mockChargeOutput()
        outputMocker.mockGammaOutput(gammaAmplitudesDict, symmetriesWeightsDict)

        # Queue slurm job
        subprocess.run(["sbatch", "job.sh"])
        os.chdir(runnerCwd)

    def createJTensorTable(self, pairingEnergySingletTripletBasis: dict[tuple[str, str], float]) -> list[float]:
        nSpins = 2
        vTensorSingletTriplet = np.zeros((nSpins * nSpins, nSpins * nSpins))
        vTensorUpDown = np.zeros((nSpins * nSpins, nSpins * nSpins))
        v4DArray = np.zeros((nSpins, nSpins, nSpins, nSpins))
        uInverse = np.array([[0,0,1,0],
                             [1/np.sqrt(2), 1/np.sqrt(2), 0, 0],
                             [-1/np.sqrt(2), 1/np.sqrt(2), 0, 0],
                             [0,0,0,1]])
        tensorSingletTripletIndexMapping: dict[str, int] = {
            "S": 0,
            "T0": 1,
            "T+": 2,
            "T-": 3,
        }
        tensorUpDownTo4DArrayMapping: dict[int, tuple[int, int]] = {
            0 : (0,0),
            1 : (0,1),
            2 : (1,0),
            3 : (1,1),
        }
        for key in pairingEnergySingletTripletBasis:
            row = tensorSingletTripletIndexMapping[key[0]]
            col = tensorSingletTripletIndexMapping[key[1]]
            if (col < row):
                raise Exception("Only upper-triangle of V-tensor should be specified")
            vTensorSingletTriplet[row][col] = pairingEnergySingletTripletBasis[key]

        indecesLowerTriangle = np.tril_indices(nSpins * nSpins, -1)
        vTensorSingletTriplet[indecesLowerTriangle] = np.conjugate(vTensorSingletTriplet.T[indecesLowerTriangle])

        vTensorUpDown = np.matmul(uInverse, np.matmul(vTensorSingletTriplet, np.transpose(uInverse)))

        for i in range(nSpins * nSpins):
            for j in range(nSpins * nSpins):
                spin1, spin2 = tensorUpDownTo4DArrayMapping[i]
                spin3, spin4 = tensorUpDownTo4DArrayMapping[j]
                v4DArray[spin1][spin2][spin3][spin4] = vTensorUpDown[i][j]

        return v4DArray.flatten(order='F').tolist()


    def __createRunDirStructure(self, runsDir: str, paramValuePairs: list[tuple[str, str, float | list[float]]]) -> str:
        pathToAppend = os.path.join(SCRATCH_PATH, runsDir)
        os.makedirs(pathToAppend, exist_ok=True)
        pathToAppend = os.path.join(pathToAppend, "RUN")

        for pair in paramValuePairs:
            if pair[0] != "self_consistency" and not isinstance(pair[2], list):
                pathToAppend = pathToAppend + f"_{pair[1]}_{pair[2]}"

        path = pathToAppend

        outputDir = f"OutputData"
        if not os.path.exists(path):
            os.mkdir(path)
            os.mkdir(os.path.join(path, outputDir))

        return path

    def __getMaterialNml(self, material: str) -> f90nml.Namelist:
        """
        Returns a namelist for self-consistent calculation that appropriate for a given material
        """
        if material == "STO":
            return self.getLaoStoDefaultNml()
        elif material == "KTO":
            return self.getLaoKtoDefaultNml()
        else:
            raise ValueError("Unknown material")

    def __createAndWriteInputNml(self, paramValuePairs: list, material: str) -> None:
        # Getting namelist with parameters
        nml = f90nml.Namelist()
        if material == "STO":
            nml = self.getLaoStoDefaultNml()
        elif material == "KTO":
            nml = self.getLaoKtoDefaultNml()

        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        with open("input.nml", "w") as nmlFile:
            f90nml.write(nml, nmlFile, sort=False)
