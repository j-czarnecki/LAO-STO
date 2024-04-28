from DataReaderClass import *
import matplotlib.pyplot as plt
import numpy as np

#TODO: self.lowestEnergy should not be used - all energies should be shown with respect to E_Fermi
#Data reader is in fact not used here, rethink this architecture
class DispersionPlotter(DataReader):

    def __init__(self):
        DataReader.__init__(self, "./", 'xxx')
        self.dataLength = 0
        self.kPoints1D = 0
        print("Initialized DispersionPlotter object")
        print(self.kPoints1D)


    def GetStatistics(self):
        self.dataLength = len(self.dispersionDataframe.N)
        self.kPoints1D = len(set(self.dispersionDataframe.kx))
        self.lowestEnergy = np.min(self.dispersionDataframe.E)

    def shiftEnergies(self):
        lowestEnergy = np.min(self.dispersionDataframe.E)
        print(lowestEnergy)
        self.dispersionDataframe.E -= lowestEnergy
        self.dosDataframe.E -= lowestEnergy

    def plotCrossection(self, plotOutputPath: str, maxEnergy: float, sliceAlong: str, fixedKVal: float, kMax: float, nBands: int):

        fixedK = ''
        xLabelOnPlot = ''
        if sliceAlong == 'kx':
            xLabelOnPlot = r'$k_x~(\tilde{a}^{-1})$'
            fixedK = 'ky'
        elif sliceAlong == 'ky':
            xLabelOnPlot = r'$k_y~(\tilde{a}^{-1})$'
            fixedK = 'kx'

        plotEnergies = np.zeros((self.kPoints1D, nBands), dtype = np.float64)
        currentIndex = np.zeros(nBands, dtype = int)
        for i in range(self.dataLength):
            if self.dispersionDataframe[fixedK][i] == fixedKVal:
                bandNo = int(self.dispersionDataframe.N[i]) - 1
                plotEnergies[currentIndex[bandNo], bandNo] = self.dispersionDataframe.E[i]
                currentIndex[bandNo] += 1
                plt.figure(0)
                plt.plot(self.dispersionDataframe[sliceAlong][i], self.dispersionDataframe.E[i], marker = '.', markersize = 2,  color = (self.dispersionDataframe.P_yz[i], self.dispersionDataframe.P_zx[i], self.dispersionDataframe.P_xy[i]))
                plt.figure(1)
                plt.plot(self.dispersionDataframe[sliceAlong][i], self.dispersionDataframe.E[i], marker = '.', markersize = 2,  color = (self.dispersionDataframe.P_lat1[i], 0, self.dispersionDataframe.P_lat2[i]))
                plt.figure(2)
                plt.plot(self.dispersionDataframe[sliceAlong][i], self.dispersionDataframe.E[i], marker = '.', markersize = 2,  color = (self.dispersionDataframe.P_up[i], 0, self.dispersionDataframe.P_down[i]))


        plt.figure(0)
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r'$E$ (meV)')
        plt.grid()
        plt.savefig(plotOutputPath + "_orbital.png")
        plt.close()

        plt.figure(1)
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r'$E$ (meV)')
        plt.grid()
        plt.savefig(plotOutputPath + "_lattice.png")
        plt.close()
       
        plt.figure(2)
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r'$E$ (meV)')
        plt.grid()
        plt.savefig(plotOutputPath + "_spin.png")
        plt.close()

        plt.figure(3)
        plt.plot(np.sort(list(set(self.dispersionDataframe.ky))), plotEnergies, linewidth = 1, color = 'black')
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r'$E$ (meV)')
        plt.grid()
        plt.savefig(plotOutputPath + '_standard.png')
        plt.close()


    def plotFermiCrossection(self, eFermi: float, dE: float, plotOutputPath: str):
        plt.figure()
        for i in range(len(self.dispersionDataframe.N)):
            if np.abs(self.dispersionDataframe.E[i] - eFermi) < dE:
                plt.plot(self.dispersionDataframe.kx[i], self.dispersionDataframe.ky[i], marker = '.', color = (self.dispersionDataframe.P_yz[i], self.dispersionDataframe.P_zx[i], self.dispersionDataframe.P_xy[i]))
        
        plt.title(r'$E_{Fermi} = $' + str(eFermi) + " (meV)")
        plt.xlabel(r'$k_x~(\tilde{a}^{-1})$')
        plt.ylabel(r'$k_y~(\tilde{a}^{-1})$')
        plt.savefig(plotOutputPath)
        plt.close()

    def plotDos(self, eMax: float, plotOutputPath: str):
        plt.figure()
        plt.plot(self.dosDataframe.DOS, self.dosDataframe.E)
        plt.ylim(bottom = 0, top=eMax)
        plt.xlabel(r'DOS')
        plt.ylabel(r"E (meV)")
        plt.grid(True)
        plt.savefig(plotOutputPath)
        plt.close()


