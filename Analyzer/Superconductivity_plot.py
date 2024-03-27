import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import f90nml


class DataReader:

    def __init__(self, runsPath, matchPattern):
        self.matchPattern = matchPattern
        self.runsPath = runsPath
        self.gamma = {}
        self.filling = {}
        self.fillingTotal = []
        self.params = []

    def __str__(self) -> str:
        dataStr = {'matchPattern': self.matchPattern,
                   'runsPath': self.runsPath}
        return str(dataStr)

    def FillDictGamma(self, pandasFile, firstIter):
        for row in range(len(pandasFile.spin)):
            dictKey = (int(pandasFile.spin[row]) , int(pandasFile.neighbor[row]) , int(pandasFile.sublat[row]) , int(pandasFile.orbital[row]))
            if firstIter:
                #self.gamma[dictKey] = [np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2)]
                self.gamma[dictKey] = [pandasFile.gammaR[row] + pandasFile.gammaIm[row]*1j]
                #print(dictKey)
            else:
                #self.gamma[dictKey].append(np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2))
                self.gamma[dictKey].append(pandasFile.gammaR[row] + pandasFile.gammaIm[row]*1j)

    def FillDictFilling(self, pandasFile, firstIter):
        self.fillingTotal.append(sum(pandasFile.filling[:]))
        for row in range(len(pandasFile.spin)):
            dictKey = (int(pandasFile.spin[row]) , int(pandasFile.sublat[row]) , int(pandasFile.orbital[row]))
            if firstIter:
                self.filling[dictKey] = [pandasFile.filling[row]]
                #print(dictKey)
            else:
                self.filling[dictKey].append(pandasFile.filling[row])
 

    def LoadFilling(self):
        directories = [dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)]
        isFirstIter = True

        for dir in directories:
            filePathConverged = os.path.join(self.runsPath, dir, 'OutputData', 'Chargen_dens_final.dat')
            filePathIter = os.path.join(self.runsPath, dir, 'OutputData', 'Chargen_dens_iter.dat')

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
        
        sortedIndexes = sorted(range(len(self.params)), key = lambda x: self.params[x])
        for key, yList in self.gamma.items():
            self.gamma[key] = [yList[i] for i in sortedIndexes]
        
        for key, yList in self.filling.items():
            self.filling[key] = [yList[i] for i in sortedIndexes]

        self.fillingTotal = [self.fillingTotal[i] for i in sortedIndexes]

        self.params = sorted(self.params)


    def LoadGamma(self, xKeywords):
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
                continue
            firstIter = False

class SymmetryResolver:
    
    def __init__(self):
        self.symmetryGamma = {}
    
    def SWavePairing(self):
        pass

    def CalculateSymmetryGamma(self, gammaDict):

        #Loop over all dict keys except from neighbours
        for spin in [1,2]:
            for sublat in [1,2]:
                for orbital in [1,2,3]:

                    #Loop over all gammas for given spin, sublat and orbital
                    for i in range(len(gammaDict[(1,1,1,1)])):
                        symmetryKey = (spin, sublat, orbital)
                        if i == 0:
                            self.symmetryGamma[symmetryKey] = self.SWavePairing(symmetryKey, gammaDict)
                        else:
                            self.symmetryGamma[symmetryKey].append(self.SWavePairing(symmetryKey, gammaDict))


def main():

    simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_U_unfinished', matchPattern= 'RUN_.*')
    simulationData.LoadFilling()
    simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    simulationData.sortData()

    #simulationDataHub = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_Hubbard', matchPattern= 'RUN_.*')
    #simulationDataHub.LoadFilling()
    #simulationDataHub.LoadGamma(xKeywords=('e_fermi', 'j_sc'))

    u_hub = [0, 2000]
    key = (1,1,1,1)
    key_filling = (1,2,1)
    #key_filling_tab = [(1,1,1), (1,2,1)]
    #plt.figure()
    for u in u_hub:
        gamma_plot = []
        n_tot_plot = []
        ef_plot = []
        n_orb_plot = []
        for i in range(len(simulationData.params)):
            if int(simulationData.params[i][1]) == int(u):
                ef_plot.append(simulationData.params[i][0])
                n_tot_plot.append(simulationData.fillingTotal[i] / 12.)
                gamma_plot.append(np.abs(simulationData.gamma[key][i]))
                n_orb_plot.append(simulationData.filling[key_filling][i])
        plt.plot(ef_plot, gamma_plot, '-', label = u)

    plt.legend(title = r"$U_{Hub}$ (meV)")
    #plt.legend()
    plt.xlabel(r"$E_{Fermi}$ (meV)")
    plt.ylabel(r"$\Gamma$ (meV)")
    plt.grid()
    #plt.xlim(0 , 0.1)
    plt.savefig("../Plots/HubbardTesting.png")



if __name__ == "__main__":
    main()