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
        self, runsPath: str, matchPattern: str, sublattices: int, subbands: int, nBands: int = 0, nAllNeighbours: int = 12
    ):
        """
        Initializes DataReader object, which contains all data from a series of simulations.
        Arguments:
            runsPath - path which contains folders with single simulations
            matchPattern - regex that tells the program which directories form runsPath should be loaded
            sublattices - number of sublattices
            subbands - number of subbands
            nBands - dimension of full Hamiltonian
        """
        # TODO: improve annotations
        self.matchPattern: str = matchPattern
        self.runsPath: str = runsPath
        self.sublattices: int = sublattices
        self.layerCouplings: int = 2 * (self.sublattices - 1)
        self.subbands: int = subbands
        self.nBands: int = nBands
        self.nAllNeighbours: int = nAllNeighbours

        # Those names are the same for k-space and real space calculations.
        self.colnamesGamma: list[str] = (
            ["band", "i_band", "j_band", "gammaR", "gammaIm"]
            #["band", "i_band", "j_band", "neighbor", "gammaR", "gammaIm"]
        )

        self.colnamesCharge: list[str] = (
            ["band", "i_band", "filling"]
        )
        self.colnamesDispersion = [
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

        self.dosColnames = ["E", "DOS", *[f"DOS_{i}" for i in range(1, self.nBands + 1)]]
        self.scGapColnames = ["kx", "ky", "gap", "state"]
        self.gammaKColnames = ["kx", "ky", "Gamma_Re", "Gamma_Im", *[f"Delta_{i}" for i in range(1, self.nAllNeighbours + 1)]]

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def LoadFilling(self, xKeywords: tuple, loadUnfinished: bool) -> pd.DataFrame:
        """
        Loads filling data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        """
        logger.info("Loading filling data")
        directories = [
            dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)
        ]

        dataframes = []

        for dir in directories:
            # TODO: it was recently corrected that filling is in Charge_dens_XXX.dat files (not Chargen_XXX)
            filePathConverged = os.path.join(
                self.runsPath, dir, "OutputData", "Charge_dens_final.dat"
            )
            filePathIter = os.path.join(
                self.runsPath, dir, "OutputData", "Charge_dens_iter.dat"
            )
            currentChargeDf = pd.DataFrame()

            if os.path.exists(filePathConverged):
                currentChargeDf = pd.read_csv(
                    filePathConverged,
                    skiprows=1,
                    comment='#',
                    sep='\s+',
                    names=self.colnamesCharge,
                    dtype=np.float64,
                )
            elif os.path.exists(filePathIter):
                logger.info(f"No convergence in {dir}")
                if loadUnfinished:
                    currentChargeDf = pd.read_csv(
                        filePathIter,
                        skiprows=1,
                        comment='#',
                        sep='\s+',
                        names=self.colnamesCharge,
                        dtype=np.float64,
                    )
            else:
                logger.info(f"No Charge dens file in {dir}")
                continue

            pathToNml = os.path.join(self.runsPath, dir)
            self.__appendXParams(pathToNml, xKeywords, currentChargeDf)
            dataframes.append(currentChargeDf)

        chargesDf = pd.concat(dataframes, ignore_index=True)
        return chargesDf

    def LoadGamma(self, xKeywords: tuple, loadUnfinished: bool) -> pd.DataFrame:
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

        dataframes = []

        for dir in directories:
            filePathGammaConverged = os.path.join(
                self.runsPath, dir, "OutputData", "Gamma_SC_final.dat"
            )
            filePathGammaIter = os.path.join(
                self.runsPath, dir, "OutputData", "Gamma_SC_iter.dat"
            )

            currentGammaDf = pd.DataFrame()
            # print(nml['physical_params'][xKeyword])
            # Gamma is printed in [meV]
            # If simulation converged final file should exists
            if os.path.exists(filePathGammaConverged):
                currentGammaDf = pd.read_csv(
                    filePathGammaConverged,
                    skiprows=1,
                    comment='#',
                    sep='\s+',
                    names=self.colnamesGamma,
                    dtype=np.float64,
                )

            # If simulation did NOT converge, iteration file should exists
            elif os.path.exists(filePathGammaIter):
                logger.info(f"No convergence in {dir}")
                if loadUnfinished:
                    currentGammaDf = pd.read_csv(
                        filePathGammaIter,
                        skiprows=1,
                        comment='#',
                        sep='\s+',
                        names=self.colnamesGamma,
                        dtype=np.float64,
                    )
            else:
                logger.warning(f"No Gamma file in {dir}")
                # shutil.rmtree(os.path.join(self.runsPath, dir))
                # print('Directory removed')
                continue

            pathToNml = os.path.join(self.runsPath, dir)
            self.__appendXParams(pathToNml, xKeywords, currentGammaDf)
            dataframes.append(currentGammaDf)

        gammasDf = pd.concat(dataframes, ignore_index=True)
        return gammasDf

    def LoadDispersion(self, energiesPath: str) -> pd.DataFrame:
        """
        Loads dispersion relations data from energiesPath.
        """
        logger.info("Loading dispersion data")
        dispersionDf = pd.DataFrame()
        if os.path.exists(energiesPath):
            dispersionDf = pd.read_csv(
                energiesPath,
                skiprows=1,
                sep='\s+',
                comment = '#',
                names=self.colnamesDispersion,
                dtype=np.float64,
            )
        else:
            logger.warning(f"No such file {energiesPath}")
        return dispersionDf

    def LoadDos(self, dosPath: str):
        """
        Loads DOS data from dosPath.
        """
        logger.info("Loading DOS data")
        dosDf = pd.DataFrame()
        if os.path.exists(dosPath):
            dosDf = pd.read_csv(
                dosPath,
                skiprows=1,
                comment='#',
                sep = '\s+',
                names=self.dosColnames,
                dtype=np.float64,
            )
        else:
            logger.warning(f"No such file {dosPath}")
        return dosDf

    def LoadSuperconductingGap(self, gapPath: str):
        """
        Loads superconducting gap from gapPath.
        """
        logger.info("Loading superconducting gap data")
        scGapDf = pd.DataFrame()
        if os.path.exists(gapPath):
            scGapDf = pd.read_csv(
                gapPath,
                skiprows=1,
                comment='#',
                sep = '\s+',
                names=self.scGapColnames,
                dtype=np.float64,
            )
        else:
            logger.warning(f"No such file {gapPath}")
        return scGapDf

    def LoadGammaMap(self, gammaKPath: str):
        gammaKDf = pd.DataFrame()
        if os.path.exists(gammaKPath):
            gammaKDf = pd.read_csv(
                gammaKPath,
                sep='\s+',
                dtype=np.float64,
                skiprows=1,
                comment='#',
                names=self.gammaKColnames,
            )
        else:
            logger.warning(f"No such file {gammaKPath}")
        return gammaKDf

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __appendXParams(self, pathToNml: str, xKeywords: tuple[str], dataframe: pd.DataFrame) -> None:
        namelistPath = os.path.join(pathToNml, "input.nml")
        with open(namelistPath) as nmlFile:
            nml = f90nml.read(nmlFile)
            for xKey in xKeywords:
                param = nml["physical_params"][xKey]
                if type(param) is list:
                    ind = [i for i, x in enumerate(param) if x != 0]
                    param = param[ind[0]]
                dataframe[xKey] = param

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __str__(self) -> str:
        dataStr = {"matchPattern": self.matchPattern, "runsPath": self.runsPath}
        return str(dataStr)
