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

import pandas as pd
import f90nml
import os
import re
import numpy as np
import shutil
import logging

logger = logging.getLogger(__name__)

class DataReader:

    def __init__(
        self, runsPath: str, matchPattern: str, sublattices: int, subbands: int
    ):
        """
        Initializes DataReader object, which contains all data from a series of simulations.
        Arguments:
            runsPath - path which contains folders with single simulations
            matchPattern - regex that tells the program which directories form runsPath should be loaded
            sublattices - number of sublattices
        Initializes:
            self.gamma - dictionary with keys of form (spin, neighbour, sublat, orbital), by default stores unsorted lists of given gammas.
                         It has the same order as self.params, based on which it could be sorted.
            self.filling - dictionary with keys of formm (spin, sublat, orbital),
                           by default contains unsorted electron concentrations in range [0,1] for each subband.
                           It has the same order as self.params, based on which it could be sorted.
            self.fillingTotal - list containing sum of all subbands' concentrations divided by 12 (number of subbands).
            self.params - contains list of tuples of parameters changed during simulations.
                          Based on that self.gamma, self.filling and self.fillingTotal can be sorted.
            self.dispersionDataFrame - contains dispersion data
            self.dosDataFrame - contains DOS data
        """
        # TODO: improve annotations
        self.matchPattern: str = matchPattern
        self.runsPath: str = runsPath
        self.sublattices: int = sublattices
        self.layerCouplings: int = 2 * (self.sublattices - 1)
        self.subbands: int = subbands
        self.gamma: dict = {}
        self.gammaKDataFrame: pd.DataFrame = pd.DataFrame()
        self.filling: dict = {}
        self.fillingTotal: list = []
        self.params: list = []
        self.dispersionDataframe = pd.DataFrame()
        self.dosDataframe = pd.DataFrame()
        self.superconductingGapDataframe = pd.DataFrame()
        self.superconductingGapMap: dict = {
            1: {"kx": [], "ky": [], "gap": []},
            2: {"kx": [], "ky": [], "gap": []},
            3: {"kx": [], "ky": [], "gap": []},
            4: {"kx": [], "ky": [], "gap": []},
            5: {"kx": [], "ky": [], "gap": []},
            6: {"kx": [], "ky": [], "gap": []},
        }
        self.colnamesGamma: list[str] = (
            ["band", "spin1", "spin2", "neighbor", "sublat", "orbital", "gammaR", "gammaIm"]
        )

        self.colnamesCharge: list[str] = (
            ["band", "spin", "sublat", "orbital", "filling"]
        )

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def LoadFilling(self, loadUnfinished: bool):
        """
        Loads filling data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        """
        logger.info("Loading filling data")
        directories = [
            dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)
        ]
        isFirstIter = True

        for dir in directories:
            # TODO: it was recently corrected that filling is in Charge_dens_XXX.dat files (not Chargen_XXX)
            filePathConverged = os.path.join(
                self.runsPath, dir, "OutputData", "Charge_dens_final.dat"
            )
            filePathIter = os.path.join(
                self.runsPath, dir, "OutputData", "Charge_dens_iter.dat"
            )
            if os.path.exists(filePathConverged):
                currentFilling = pd.read_fwf(
                    filePathConverged,
                    skiprows=1,
                    names=self.colnamesCharge,
                    infer_nrows=10,
                    dtype=np.float64,
                )
                self.__FillDictFilling(currentFilling, isFirstIter)
            elif os.path.exists(filePathIter):
                logger.info(f"No convergence in {dir}")
                if loadUnfinished:
                    currentFilling = pd.read_fwf(
                        filePathIter,
                        skiprows=1,
                        names=self.colnamesCharge,
                        infer_nrows=10,
                        dtype=np.float64,
                    )
                    self.__FillDictFilling(currentFilling, isFirstIter)
                else:
                    self.__FillDictFilling(currentFilling, isFirstIter, fillNones=True)
            else:
                logger.info(f"No Charge dens file in {dir}")
                continue
            isFirstIter = False

    def LoadGamma(self, xKeywords: tuple, loadUnfinished: bool):
        """
        Loads gamma data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        Additionally fills self.params list based on xKeywords - names of f90 .nml parameters from input.nml file
        that were changed during simulation.
        """

        logger.info("Loading gamma data")
        directories = [
            dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)
        ]

        firstIter = True

        for dir in directories:
            filePathGammaConverged = os.path.join(
                self.runsPath, dir, "OutputData", "Gamma_SC_final.dat"
            )
            filePathGammaIter = os.path.join(
                self.runsPath, dir, "OutputData", "Gamma_SC_iter.dat"
            )

            # print(nml['physical_params'][xKeyword])
            # Gamma is printed in [meV]
            # If simulation converged final file should exists
            if os.path.exists(filePathGammaConverged):
                currentGamma = pd.read_fwf(
                    filePathGammaConverged,
                    skiprows=1,
                    names=self.colnamesGamma,
                    infer_nrows=100,
                    dtype=np.float64,
                )
                self.__FillDictGamma(currentGamma, firstIter)

            # If simulation did NOT converge, iteration file should exists
            elif os.path.exists(filePathGammaIter):
                logger.info(f"No convergence in {dir}")
                if loadUnfinished:
                    currentGamma = pd.read_fwf(
                        filePathGammaIter,
                        skiprows=1,
                        names=self.colnamesGamma,
                        infer_nrows=100,
                        dtype=np.float64,
                    )
                    self.__FillDictGamma(currentGamma, firstIter)
                else:
                    self.__FillDictGamma(currentGamma, firstIter, fillNones=True)
            else:
                logger.warning(f"No Gamma file in {dir}")
                # shutil.rmtree(os.path.join(self.runsPath, dir))
                # print('Directory removed')
                continue

            namelistPath = os.path.join(self.runsPath, dir, "input.nml")
            with open(namelistPath) as nmlFile:
                nml = f90nml.read(nmlFile)
                paramsValuesList = []
                for xKey in xKeywords:
                    param = nml["physical_params"][xKey]
                    if type(param) is list:
                        ind = [i for i, x in enumerate(param) if x != 0]
                        param = param[ind[0]]
                    paramsValuesList.append(param)
                self.params.append(tuple(paramsValuesList))

            firstIter = False

    def LoadDispersion(self, energiesPath: str):
        """
        Loads dispersion relations data from energiesPath.
        """
        logger.info("Loading dispersion data")
        names = [
            "N",
            "kx",
            "ky",
            "E",
            "P_yz",
            "P_zx",
            "P_xy",
            *[f"P_lat{i}" for i in range(1, self.sublattices + 1)],
            "P_sx",
            "P_sy",
            "P_sz",
            "P_elec",
            "P_hole",
        ]
        real_width = 14
        cols = [
            (0, 6),
            *[
                (7 + i + i * real_width, 7 + i + (i + 1) * real_width)
                for i in range(11 + self.sublattices)
            ],
        ]
        if os.path.exists(energiesPath):
            self.dispersionDataframe = pd.read_fwf(
                energiesPath,
                skiprows=1,
                colspecs=cols,
                names=names,
                dtype=np.float32,
                low_memory=True,
            )
        else:
            logger.warning(f"No such file {energiesPath}")

    def LoadDos(self, dosPath: str):
        """
        Loads DOS data from dosPath.
        """
        logger.info("Loading DOS data")
        if os.path.exists(dosPath):
            self.dosDataframe = pd.read_fwf(
                dosPath,
                skiprows=1,
                colspecs=[(0, 16), (17, 31)],
                names=["E", "DOS"],
                dtype=np.float64,
            )
        else:
            logger.warning(f"No such file {dosPath}")

    def LoadSuperconductingGap(self, gapPath: str):
        """
        Loads superconducting gap from gapPath.
        """
        logger.info("Loading superconducting gap data")
        if os.path.exists(gapPath):
            self.superconductingGapDataframe = pd.read_fwf(
                gapPath,
                skiprows=1,
                infer_nrows=-100,
                colspecs=[(0, 16), (17, 31), (32, 46), (47, 56)],
                names=["kx", "ky", "gap", "state"],
                dtype=np.float64,
            )
        else:
            logger.warning(f"No such file {gapPath}")

    def LoadSuperconductingGapMap(self, runsPathGap: str, matchPatternGap: str):
        directories = [
            dir for dir in os.listdir(runsPathGap) if re.match(matchPatternGap, dir)
        ]

        for dir in directories:
            filePath = os.path.join(
                runsPathGap, dir, "OutputData", "SuperconductingGap.dat"
            )
            if os.path.exists(filePath):
                currentGap = pd.read_fwf(
                    filePath,
                    skiprows=1,
                    infer_nrows=-100,
                    colspecs=[(0, 16), (17, 31), (32, 46), (47, 56)],
                    names=["kx", "ky", "gap", "state"],
                    dtype=np.float64,
                )
                self.__fillDictScGap(currentGap)
            else:
                logger.warning(f"No Gap file in {dir}")

    def LoadGammaMap(self, gammaKPath: str):
        if os.path.exists(gammaKPath):
            self.gammaKDataFrame = pd.read_csv(
                gammaKPath,
                delim_whitespace=True,
                dtype=np.float64,
                skiprows=1,
            )
        else:
            logger.warning(f"No such file {gammaKPath}")

    def sortData(self):
        """
        Sorts datafrom self.gamma, self.filling and self.fillingTotal dicts.
        Based on the order of self.params list of tuples (first elements ???)
        """
        sortedIndexes = sorted(range(len(self.params)), key=lambda x: self.params[x])
        for key, yList in self.gamma.items():
            self.gamma[key] = [yList[i] for i in sortedIndexes]

        for key, yList in self.filling.items():
            self.filling[key] = [yList[i] for i in sortedIndexes]

        self.fillingTotal = [self.fillingTotal[i] for i in sortedIndexes]

        self.params = sorted(self.params)

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __FillDictGamma(
        self, pandasFile: pd.DataFrame, firstIter: bool, fillNones: bool = False
    ):
        """
        Extracts proper key from simulation files, converts gamma back to complex value and writes do dict
        """
        if len(pandasFile.spin1) == 0:
            logger.error("Empty gamma file")
            for key in list(self.gamma.keys()):
                self.gamma[key].append(np.nan)
            return

        for row in range(len(pandasFile["spin1"])):
            # Key is everything apart from Gamma values (two last columns)
            dictKey = tuple(
                int(x) for x in pandasFile.loc[row, self.colnamesGamma[:-2]]
            )
            if firstIter:
                if not fillNones:
                    self.gamma[dictKey] = [
                        pandasFile.gammaR[row] + pandasFile.gammaIm[row] * 1j
                    ]
                else:
                    self.gamma[dictKey] = np.nan
            else:
                if not fillNones:
                    self.gamma[dictKey].append(
                        pandasFile.gammaR[row] + pandasFile.gammaIm[row] * 1j
                    )
                else:
                    self.gamma[dictKey].append(np.nan)

    def __FillDictFilling(
        self, pandasFile: pd.DataFrame, firstIter: bool, fillNones: bool = False
    ):
        """
        Extracts proper key from simulation files and appends to filling dict
        """
        if len(pandasFile.spin) == 0:
            logger.error("Empty filling file")
            for key in list(self.filling.keys()):
                self.filling[key].append(np.nan)
            self.fillingTotal.append(np.nan)
            return

        self.fillingTotal.append(sum(pandasFile.filling[:]))
        for row in range(len(pandasFile.spin)):
            dictKey = tuple(
                int(x) for x in pandasFile.loc[row, self.colnamesCharge[:-1]]
            )

            if firstIter:
                if not fillNones:
                    self.filling[dictKey] = [pandasFile.filling[row]]
                else:
                    self.filling[dictKey] = np.nan
            else:
                if not fillNones:
                    self.filling[dictKey].append(pandasFile.filling[row])
                else:
                    self.filling[dictKey].append(np.nan)

    def __fillDictScGap(self, pandasFile: pd.DataFrame, fillNones: bool = False):
        for row in range(len(pandasFile["kx"])):
            self.superconductingGapMap[int(pandasFile["state"][row])]["kx"].append(
                pandasFile["kx"][row]
            )
            self.superconductingGapMap[int(pandasFile["state"][row])]["ky"].append(
                pandasFile["ky"][row]
            )
            self.superconductingGapMap[int(pandasFile["state"][row])]["gap"].append(
                pandasFile["gap"][row]
            )

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __str__(self) -> str:
        dataStr = {"matchPattern": self.matchPattern, "runsPath": self.runsPath}
        return str(dataStr)
