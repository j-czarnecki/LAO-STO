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



#This class should be improved
class SymmetryResolver:
    
    def __init__(self, nNeighbors):
        self.nNeighbors = nNeighbors
        self.symmetryGammaDict = {}

    #RETHINK ___WavePairing construction to avoid code repetition
    def _SWavePairing(self, listOfGammas: list):
        p = 0
        M = 0
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _PWavePairing(self, listOfGammas: list):
        p = 1
        M = 1
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _DWavePairing(self, listOfGammas: list):
        p = 0
        M = 2
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _CalculateSingleGamma(self, listOfGammas: list, p: int, M: int):
        symmetryGamma = 0
        for i in range(self.nNeighbors):
            currentPhase = 2*np.pi / self.nNeighbors * i
            phaseFactor = np.exp(-1j*M*currentPhase)
            symmetryGamma += listOfGammas[i]*phaseFactor

        return symmetryGamma*(1j)**p/self.nNeighbors

    def CalculateSymmetryGamma(self, gammaDict):

        #Loop over all dict keys except from neighbours
        for spin in [1,2]:
            for sublat in [1,2]:
                for orbital in [1,2,3]:

                    #Loop over all gammas for given spin, sublat and orbital
                    for i in range(len(gammaDict[(1,1,1,1)])):
                        gammaToSymmetrize = []
                        for neighbor in range(1, self.nNeighbors + 1):
                            gammaToSymmetrize.append(gammaDict[(spin, neighbor, sublat, orbital)][i])
                        if i == 0:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')] = [self._SWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p')] = [self._PWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd')] = [self._DWavePairing(gammaToSymmetrize)]
                        else:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')].append(self._SWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p')].append(self._PWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd')].append(self._DWavePairing(gammaToSymmetrize))

def testing():
    symmetryResolver = SymmetryResolver(3)

    testGammaDict = {}
    for spin in [1,2]:
        for sublat in [1,2]:
            for orbital in [1,2,3]:
                for neighbor in [1,2,3]:
                    key = (spin, neighbor, sublat, orbital)
                    testGammaDict[key] = [x*neighbor for x in range(1,3)]

    print(testGammaDict)
    symmetryResolver.CalculateSymmetryGamma(testGammaDict)
    print(symmetryResolver.symmetryGammaDict[(1,1,1,'s')])
    print(symmetryResolver.symmetryGammaDict[(1,1,1,'p')])
    print(symmetryResolver.symmetryGammaDict[(1,1,1,'d')])

def plotGammasFermi(simulationData, symmetryResolver):
    U_tab = [0, 166, 333]
    j_sc_tab = [75,150]
    keys = []

    for spin in [1]: #should be [1,2], but now spin are symmetric
        for sublat in [1]: #should be [1,2], but now sublattices are symmetric
            for orbital in [1,2,3]:
                for symmetry in ['s','p','d']:
                    keys.append((spin, sublat, orbital, symmetry))


    for key in keys:
        plt.figure()
        for u in U_tab:
            gamma_plot = []
            ef_plot= []

            for i in range(len(simulationData.params)):
                if int(simulationData.params[i][1]) == u:
                    ef_plot.append(simulationData.params[i][0])
                    gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
            plt.plot(ef_plot, gamma_plot, '-', label = u)

        spin, sublat, orbital, symmetry = key
        plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
        plt.legend(title = r"$U_{Hub}$ (meV)")
        plt.xlabel(r"$E_{Fermi}$ (meV)")
        plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
        plt.grid()
        #plt.xlim(0 , 0.1)
        plt.savefig(f"../Plots/GammaFermi_{spin}_{sublat}_{orbital}_{symmetry}.png")
        plt.close()

def plotGammasFilling(simulationData, symmetryResolver):
    U_tab = [0, 166, 333]
    j_sc_tab = [75,150]
    keys = []

    for spin in [1]: #should be [1,2], but now spin are symmetric
        for sublat in [1]: #should be [1,2], but now sublattices are symmetric
            for orbital in [1,2,3]:
                for symmetry in ['s','p','d']:
                    keys.append((spin, sublat, orbital, symmetry))


    for key in keys:
        plt.figure()
        for u in U_tab:
            gamma_plot = []
            n_total_plot = []

            for i in range(len(simulationData.params)):
                if int(simulationData.params[i][1]) == u:
                    n_total_plot.append(simulationData.fillingTotal[i]/12.)
                    #gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                    gamma_plot.append(np.abs(simulationData.gamma[(2,2,2,2)][i]))
            plt.plot(n_total_plot, gamma_plot, '-', label = u)

        spin, sublat, orbital, symmetry = key
        plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
        plt.legend(title = r"$U_{Hub}$ (meV)")
        plt.xlabel(r"$n_{tot}$")
        plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
        plt.grid()
        #plt.xlim(0 , 0.1)
        plt.savefig(f"../Plots/GammaFilling_{spin}_{sublat}_{orbital}_{symmetry}.png")
        plt.close() 

def plotFillingFermi(simulationData, symmetryResolver):
    U_tab = [0, 166, 333]
    plt.figure()
    for u in U_tab:
        ef_plot = []
        n_total_plot = []
        n_chosen_plot = []

        for i in range(len(simulationData.params)):
            if int(simulationData.params[i][1]) == u:
                n_total_plot.append(simulationData.fillingTotal[i]/12.)
                n_chosen_plot.append(simulationData.filling[(1,2,1)][i]/12.) #key (spin,sublat,orbital)
                #gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                ef_plot.append(simulationData.params[i][0])
        plt.plot(ef_plot, n_chosen_plot, '-', label = u)

    plt.legend(title = r"$U_{Hub}$ (meV)")
    plt.xlabel(r"$E_{Fermi}$")
    plt.ylabel(r"$n_{tot}$")
    plt.grid()
    #plt.xlim(0 , 0.1)
    plt.savefig(f"../Plots/FillingFermi.png")
    plt.close() 

def main():

    simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_low_U', matchPattern= 'RUN_.*')
    simulationData.LoadFilling()
    simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    simulationData.sortData()

    symmetryResolver = SymmetryResolver(3)
    symmetryResolver.CalculateSymmetryGamma(simulationData.gamma)

    plotGammasFermi(simulationData, symmetryResolver)
    plotGammasFilling(simulationData, symmetryResolver)
    plotFillingFermi(simulationData, symmetryResolver)

    # simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_U_unfinished', matchPattern= 'RUN_.*')
    # simulationData.LoadFilling()
    # simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    # simulationData.sortData()

    #simulationDataHub = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_Hubbard', matchPattern= 'RUN_.*')
    #simulationDataHub.LoadFilling()
    #simulationDataHub.LoadGamma(xKeywords=('e_fermi', 'j_sc'))



if __name__ == "__main__":
    main()
    #testing()