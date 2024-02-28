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
                self.gamma[dictKey] = [np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2)]
                #print(dictKey)
            else:
                self.gamma[dictKey].append(np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2))

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
                
        self.fillingTotal = sorted(self.fillingTotal)

    def sortDict(self):
        
        sortedIndexes = sorted(range(len(self.params)), key = lambda x: self.params[x])
        for key, yList in self.gamma.items():
            self.gamma[key] = [yList[i] for i in sortedIndexes]
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

        self.sortDict()

def main():

    simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_EF_refined', matchPattern= 'RUN_.*')
    simulationData.LoadFilling()
    simulationData.LoadGamma(xKeywords=('e_fermi', 'j_sc'))


    simulationDataHub = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_Hubbard', matchPattern= 'RUN_.*')
    simulationDataHub.LoadFilling()
    simulationDataHub.LoadGamma(xKeywords=('e_fermi', 'j_sc'))

    j_sc = 150.
    key = (1,1,1,1)


    #Comparison of hubbard and without hubbard simulation
    X_no_hub = []
    Y_no_hub = []
    X_hub = []
    Y_hub = []
    plt.figure()
    for i in range(len(simulationData.params)):
        if int(simulationData.params[i][1]) == int(j_sc):
            X_no_hub.append(simulationData.params[i][0])
            #X_no_hub.append(simulationData.fillingTotal[i])
            Y_no_hub.append(simulationData.gamma[key][i])
    
    for i in range(len(simulationDataHub.params)):
        if int(simulationDataHub.params[i][1]) == int(j_sc):
            X_hub.append(simulationDataHub.params[i][0])
            #X_hub.append(simulationDataHub.fillingTotal[i])
            Y_hub.append(simulationDataHub.gamma[key][i])
    

    plt.title(r"$J_{SC} = 150$ (meV)")
    plt.plot(X_no_hub, Y_no_hub, '-', label = ' U = V = 0')
    plt.plot(X_hub, Y_hub, '-', label = 'U = V = 2000 meV')
    plt.xlabel(r"$E_{Fermi}$ (meV)")
    plt.ylabel(r"$\Gamma$ (meV)")
    plt.legend()
    plt.savefig("../Plots/GammaHubComp.png")


    #As a function of E_fermi for different J_SC
    j_tab = [75.25, 150.]
    plt.figure()
    for j_sc in j_tab:
        X_plot = []
        Y_plot = []        
        for i in range(len(simulationData.params)):
            if int(simulationData.params[i][1]) == int(j_sc):
                X_plot.append(simulationData.params[i][0])
                #X_hub.append(simulationData.fillingTotal[i])
                Y_plot.append(simulationData.gamma[key][i])
        plt.plot(X_plot, Y_plot, '-', label = j_sc)
    
    plt.title(r"$U = V = 0$")
    plt.xlabel(r"$E_{Fermi}$ (meV)")
    plt.ylabel(r"$\Gamma$ (meV)")
    plt.legend(title = r"$J_{SC}$ (meV)")
    plt.savefig("../Plots/GammaRefined.png")

    



if __name__ == "__main__":
    main()