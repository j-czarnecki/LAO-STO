import pandas as pd
import f90nml
import os
import re
import numpy as np
import shutil

class DataReader:

    def __init__(self, runsPath: str, matchPattern: str):
        '''
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
        '''
        #TODO: improve annotations
        self.matchPattern = matchPattern
        self.runsPath = runsPath
        self.gamma: dict = {}
        self.filling: dict = {}
        self.fillingTotal: list = []
        self.params: list = []
        self.dispersionDataframe = pd.DataFrame()
        self.dosDataframe = pd.DataFrame()

    def __str__(self) -> str:
        dataStr = {'matchPattern': self.matchPattern,
                   'runsPath': self.runsPath}
        return str(dataStr)

    def FillDictGamma(self, pandasFile: pd.DataFrame, firstIter: bool):
        """
        Extracts proper key from simulation files, converts gamma back to complex value and writes do dict
        """
        for row in range(len(pandasFile.spin)):
            dictKey = (int(pandasFile.spin[row]) , int(pandasFile.neighbor[row]) , int(pandasFile.sublat[row]) , int(pandasFile.orbital[row]))
            if firstIter:
                #self.gamma[dictKey] = [np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2)]
                self.gamma[dictKey] = [pandasFile.gammaR[row] + pandasFile.gammaIm[row]*1j]
                #print(dictKey)
            else:
                #self.gamma[dictKey].append(np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2))
                self.gamma[dictKey].append(pandasFile.gammaR[row] + pandasFile.gammaIm[row]*1j)

    def FillDictFilling(self, pandasFile: pd.DataFrame, firstIter: bool):
        """
        Extracts proper key from simulation files and appends to filling dict
        """
        self.fillingTotal.append(sum(pandasFile.filling[:]))
        for row in range(len(pandasFile.spin)):
            dictKey = (int(pandasFile.spin[row]) , int(pandasFile.sublat[row]) , int(pandasFile.orbital[row]))
            if firstIter:
                self.filling[dictKey] = [pandasFile.filling[row]]
                #print(dictKey)
            else:
                self.filling[dictKey].append(pandasFile.filling[row])
 

    def LoadFilling(self):
        """
        Loads filling data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        """
        print("---> Loading filling data")
        directories = [dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)]
        isFirstIter = True

        for dir in directories:
            #TODO: it was recently corrected that filling is in Charge_dens_XXX.dat files (not Chargen_XXX)
            filePathConverged = os.path.join(self.runsPath, dir, 'OutputData', 'Charge_dens_final.dat')
            filePathIter = os.path.join(self.runsPath, dir, 'OutputData', 'Charge_dens_iter.dat')

            if os.path.exists(filePathConverged):
                currentFilling = pd.read_fwf(filePathConverged, skiprows = 1, colspecs = [(0,6), (7,11), (12,16), (17,31)],
                                             names = ['spin', 'sublat', 'orbital', 'filling'])
                self.FillDictFilling(currentFilling, isFirstIter)
            elif os.path.exists(filePathIter):
                currentFilling = pd.read_fwf(filePathIter, skiprows = 1, colspecs = [(0,6), (7,11), (12,16), (17,31)],
                                             names = ['spin', 'sublat', 'orbital', 'filling'])
                self.FillDictFilling(currentFilling, isFirstIter)
            else:
                print("No Charge dens file in ", dir)
                continue
            isFirstIter = False
                
    def sortData(self):
        """
        Sorts datafrom self.gamma, self.filling and self.fillingTotal dicts.
        Based on the order of self.params list of tuples (first elements ???)
        """
        sortedIndexes = sorted(range(len(self.params)), key = lambda x: self.params[x])
        for key, yList in self.gamma.items():
            self.gamma[key] = [yList[i] for i in sortedIndexes]
        
        for key, yList in self.filling.items():
            self.filling[key] = [yList[i] for i in sortedIndexes]

        self.fillingTotal = [self.fillingTotal[i] for i in sortedIndexes]

        self.params = sorted(self.params)


    def LoadGamma(self, xKeywords: tuple):
        """
        Loads gamma data from simulations base on specified in __init__() runsPath and matchPattern.
        If simulation had not converged, takse values from _iter.dat file - the last iteration before program timeout.
        Additionally fills self.params list based on xKeywords - names of f90 .nml parameters from input.nml file
        that were changed during simulation.
        """

        print("---> Loading gamma data")
        directories = [dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)]
        #print(directories)

        firstIter = True

        for dir in directories:
            filePathGammaConverged = os.path.join(self.runsPath, dir, 'OutputData', 'Gamma_SC_final.dat')
            filePathGammaIter = os.path.join(self.runsPath, dir, 'OutputData', 'Gamma_SC_iter.dat')
            namelistPath = os.path.join(self.runsPath, dir, 'input.nml')
            
            with open(namelistPath) as nmlFile:
                nml = f90nml.read(nmlFile)
                paramsValuesList = []
                for xKey in xKeywords:
                    paramsValuesList.append(nml['physical_params'][xKey])
                self.params.append(tuple(paramsValuesList))
                #print(nml['physical_params'][xKeyword])
            
            #Gamma is printed in [meV]
            #If simulation converged final file should exists
            if os.path.exists(filePathGammaConverged):        
                currentGamma = pd.read_fwf(filePathGammaConverged, skiprows = 1, colspecs = [(0,6), (7,11), (12, 16), (17,21), (22,36), (37,51)], \
                                        names = ['spin', 'neighbor', 'sublat', 'orbital', 'gammaR', 'gammaIm'], dtype = np.float64)
                self.FillDictGamma(currentGamma, firstIter)   
                
            #If simulation did NOT converge, iteration file should exists
            elif os.path.exists(filePathGammaIter):        
                currentGamma = pd.read_fwf(filePathGammaIter, skiprows = 1, colspecs = [(0,6), (7,11), (12, 16), (17,21), (22,36), (37,51)], names = ['spin', 'neighbor', 'sublat', 'orbital', 'gammaR', 'gammaIm'], dtype = np.float64)
                self.FillDictGamma(currentGamma, firstIter)
                print("No convergence in ", dir)       
            
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
           #dispersion = pd.read_fwf(energiesPath, skiprows=1, colspecs=[(0,6), (7,21), (22,36), (37,51), (52, 66), (67, 81), (82, 96), (97, 111), (112,126), (127, 141), (142, 156)]) 
            self.dispersionDataframe = pd.read_fwf(energiesPath, skiprows = 1, infer_nrows=-100,
                                     colspecs=[(0,6), (7,21), (22,36), (37,51), (52, 66), (67, 81), (82, 96), (97, 111), (112,126), (127, 141), (142, 156)],
                                     names = ['N', 'kx', 'ky', 'E', 'P_yz', 'P_zx', 'P_xy', 'P_lat1', 'P_lat2', 'P_up', 'P_down'], dtype = np.float64 )
        else:
            print("No such file ", energiesPath)
        

    def LoadDos(self, dosPath: str):
        """
        Loads DOS data from dosPath.
        """
        print("---> Loading DOS data")
        if os.path.exists(dosPath):
           #dispersion = pd.read_fwf(energiesPath, skiprows=1, colspecs=[(0,6), (7,21), (22,36), (37,51), (52, 66), (67, 81), (82, 96), (97, 111), (112,126), (127, 141), (142, 156)]) 
            self.dosDataframe = pd.read_fwf(dosPath, skiprows = 1, infer_nrows=-100, colspecs=[(0,16), (17,31)], names = ['E', 'DOS'], dtype = np.float64 )
        else:
            print("No such file ", dosPath)

