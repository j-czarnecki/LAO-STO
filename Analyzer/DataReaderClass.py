import pandas as pd
import f90nml
import os
import re
import numpy as np
import shutil


class DataReader:

    def __init__(self, runsPath: str, matchPattern: str):
        """
        Initializes DataReader object, which contains all data from a series of simulations.
        Arguments:
            runsPath - path which contains folders with single simulations
            matchPattern - regex that tells the program which directories form runsPath should be loaded
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
        self.matchPattern = matchPattern
        self.runsPath = runsPath
        self.gamma: dict = {}
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

    def __str__(self) -> str:
        dataStr = {"matchPattern": self.matchPattern, "runsPath": self.runsPath}
        return str(dataStr)

    def FillDictGamma(
        self, pandasFile: pd.DataFrame, firstIter: bool, fillNones: bool = False
    ):
        """
        Extracts proper key from simulation files, converts gamma back to complex value and writes do dict
        """
        if len(pandasFile.spin) == 0:
            print("ERROR: empty gamma file")
            for key in list(self.gamma.keys()):
                self.gamma[key].append(np.nan)
            return

        for row in range(len(pandasFile.spin)):
            dictKey = (
                int(pandasFile.spin[row]),
                int(pandasFile.neighbor[row]),
                int(pandasFile.sublat[row]),
                int(pandasFile.orbital[row]),
            )
            if firstIter:
                # self.gamma[dictKey] = [np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2)]
                if not fillNones:
                    self.gamma[dictKey] = [
                        pandasFile.gammaR[row] + pandasFile.gammaIm[row] * 1j
                    ]
                else:
                    self.gamma[dictKey] = np.nan
                # print(dictKey)
            else:
                # self.gamma[dictKey].append(np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2))
                if not fillNones:
                    self.gamma[dictKey].append(
                        pandasFile.gammaR[row] + pandasFile.gammaIm[row] * 1j
                    )
                else:
                    self.gamma[dictKey].append(np.nan)

    def FillDictFilling(
        self, pandasFile: pd.DataFrame, firstIter: bool, fillNones: bool = False
    ):
        """
        Extracts proper key from simulation files and appends to filling dict
        """
        if len(pandasFile.spin) == 0:
            print("ERROR: empty filling file")
            for key in list(self.filling.keys()):
                self.filling[key].append(np.nan)
            self.fillingTotal.append(np.nan)
            return

        self.fillingTotal.append(sum(pandasFile.filling[:]))
        for row in range(len(pandasFile.spin)):
            dictKey = (
                int(pandasFile.spin[row]),
                int(pandasFile.sublat[row]),
                int(pandasFile.orbital[row]),
            )
            if firstIter:
                if not fillNones:
                    self.filling[dictKey] = [pandasFile.filling[row]]
                else:
                    self.filling[dictKey] = np.nan
                # print(dictKey)
            else:
                if not fillNones:
                    self.filling[dictKey].append(pandasFile.filling[row])
                else:
                    self.filling[dictKey].append(np.nan)

    def FillDictScGap(self, pandasFile: pd.DataFrame, fillNones: bool = False):
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

    def LoadFilling(self, loadUnfinished: bool):
        """
        Loads filling data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        """
        print("---> Loading filling data")
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
                    colspecs=[(0, 6), (7, 11), (12, 16), (17, 31)],
                    names=["spin", "sublat", "orbital", "filling"],
                )
                self.FillDictFilling(currentFilling, isFirstIter)
            elif os.path.exists(filePathIter):
                print("No convergence in ", dir)
                if loadUnfinished:
                    currentFilling = pd.read_fwf(
                        filePathIter,
                        skiprows=1,
                        colspecs=[(0, 6), (7, 11), (12, 16), (17, 31)],
                        names=["spin", "sublat", "orbital", "filling"],
                    )
                    self.FillDictFilling(currentFilling, isFirstIter)
                else:
                    self.FillDictFilling(currentFilling, isFirstIter, fillNones=True)
            else:
                print("No Charge dens file in ", dir)
                continue
            isFirstIter = False

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

    def LoadGamma(self, xKeywords: tuple, loadUnfinished: bool):
        """
        Loads gamma data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        Additionally fills self.params list based on xKeywords - names of f90 .nml parameters from input.nml file
        that were changed during simulation.
        """

        print("---> Loading gamma data")
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
            namelistPath = os.path.join(self.runsPath, dir, "input.nml")

            with open(namelistPath) as nmlFile:
                nml = f90nml.read(nmlFile)
                paramsValuesList = []
                for xKey in xKeywords:
                    paramsValuesList.append(nml["physical_params"][xKey])
                self.params.append(tuple(paramsValuesList))
                # print(nml['physical_params'][xKeyword])

            # Gamma is printed in [meV]
            # If simulation converged final file should exists
            if os.path.exists(filePathGammaConverged):
                currentGamma = pd.read_fwf(
                    filePathGammaConverged,
                    skiprows=1,
                    colspecs=[(0, 6), (7, 11), (12, 16), (17, 21), (22, 36), (37, 51)],
                    names=[
                        "spin",
                        "neighbor",
                        "sublat",
                        "orbital",
                        "gammaR",
                        "gammaIm",
                    ],
                    dtype=np.float64,
                )
                self.FillDictGamma(currentGamma, firstIter)

            # If simulation did NOT converge, iteration file should exists
            elif os.path.exists(filePathGammaIter):
                print("No convergence in ", dir)
                if loadUnfinished:
                    currentGamma = pd.read_fwf(
                        filePathGammaIter,
                        skiprows=1,
                        colspecs=[
                            (0, 6),
                            (7, 11),
                            (12, 16),
                            (17, 21),
                            (22, 36),
                            (37, 51),
                        ],
                        names=[
                            "spin",
                            "neighbor",
                            "sublat",
                            "orbital",
                            "gammaR",
                            "gammaIm",
                        ],
                        dtype=np.float64,
                    )
                    self.FillDictGamma(currentGamma, firstIter)
                else:
                    self.FillDictGamma(currentGamma, firstIter, fillNones=True)
            else:
                print("No Gamma file in ", dir)
                # shutil.rmtree(os.path.join(self.runsPath, dir))
                # print('Direcotry removed')
                continue
            firstIter = False

    def LoadDispersion(self, energiesPath: str):
        """
        Loads dispersion relations data from energiesPath.
        """
        print("---> Loading dispersion data")
        if os.path.exists(energiesPath):
            self.dispersionDataframe = pd.read_fwf(
                energiesPath,
                skiprows=1,
                colspecs=[
                    (0, 6),
                    (7, 21),
                    (22, 36),
                    (37, 51),
                    (52, 66),
                    (67, 81),
                    (82, 96),
                    (97, 111),
                    (112, 126),
                    (127, 141),
                    (142, 156),
                    (157, 171),
                    (172, 186),
                ],
                names=[
                    "N",
                    "kx",
                    "ky",
                    "E",
                    "P_yz",
                    "P_zx",
                    "P_xy",
                    "P_lat1",
                    "P_lat2",
                    "P_up",
                    "P_down",
                    "P_elec",
                    "P_hole",
                ],
                dtype=np.float32,
                low_memory=True,
            )
        else:
            print("No such file ", energiesPath)


    def LoadDos(self, dosPath: str):
        """
        Loads DOS data from dosPath.
        """
        print("---> Loading DOS data")
        if os.path.exists(dosPath):
            # dispersion = pd.read_fwf(energiesPath, skiprows=1, colspecs=[(0,6), (7,21), (22,36), (37,51), (52, 66), (67, 81), (82, 96), (97, 111), (112,126), (127, 141), (142, 156)])
            self.dosDataframe = pd.read_fwf(
                dosPath,
                skiprows=1,
                colspecs=[(0, 16), (17, 31)],
                names=["E", "DOS"],
                dtype=np.float64,
            )
        else:
            print("No such file ", dosPath)

    def LoadSuperconductingGap(self, gapPath: str):
        """
        Loads superconducting gap from gapPath.
        """
        print("---> Loading superconducting gap")
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
            print("No such file ", gapPath)

    def LoadSuperconductingGapMap(self, runsPathGap: str, matchPatternGap: str):
        directories = [
            dir for dir in os.listdir(runsPathGap) if re.match(matchPatternGap, dir)
        ]
        isFirstIter = True

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
                self.FillDictScGap(currentGap)
            else:
                print("No such file ", filePath)
            isFirstIter = False
