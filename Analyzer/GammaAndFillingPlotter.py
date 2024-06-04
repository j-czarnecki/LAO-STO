import matplotlib.pyplot as plt
import numpy as np
from DataReaderClass import *
from SymmetryResolverClass import *

#TODO: this class should be improved to be more general and possibly plot more symmetries of gamma
# self.eMinimal should not be used, as all energies must be calculated with respect to E_Fermi
class GammaAndFillingPlotter(SymmetryResolver):

    def __init__(self, runsPath: str, matchPattern: str, nNeighbors: int, nNextNeighbors: int, eMinimal: float):
        SymmetryResolver.__init__(self, nNeighbors, nNextNeighbors, runsPath, matchPattern)
        self.eMinimal = eMinimal
        self.symmetryKeys: list[tuple[int,int,int,str]] = []
        self.nnnSymmetryKeys: list[tuple[int,int,int,str]] = []
        self.__initializedSymmetryKeys()
        print("Initialized GammaAndFillingPlotter")

    def __initializedSymmetryKeys(self):
        for spin in (1,2): #should be [1,2], but now spin are symmetric
            for sublat in (1,2): #should be [1,2], but now sublattices are symmetric
                for orbital in (1,2,3):
                    for symmetry in ('s','p+','d+'):
                        self.symmetryKeys.append((spin, sublat, orbital, symmetry))

        for spin in (1,2): #should be [1,2], but now spin are symmetric
            for sublat in (1,2): #should be [1,2], but now sublattices are symmetric
                for orbital in (1,2,3):
                    for symmetry in ('s','p+','d+', 'f', 'd-', 'p-'):
                        self.nnnSymmetryKeys.append((spin, sublat, orbital, symmetry))

    
    def plotGammasJ(self):
        E_fermi_tab = [-1030, -1040, -1050]

        for key in self.symmetryKeys:
            plt.figure()
            for ef in E_fermi_tab:
                gamma_plot = []
                j_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == ef:
                        j_plot.append(self.params[i][0] - self.eMinimal)
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                plt.plot(j_plot, gamma_plot, '-', label = ef)

            spin, sublat, orbital, symmetry = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$E_{Fermi}$ (meV)")
            plt.xlabel(r"$J_{SC}$ (meV)")
            plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
            plt.grid()
            plt.ylim(0 , 0.2)
            plt.savefig(f"../Plots/GammaJ_{spin}_{sublat}_{orbital}_{symmetry}.png")
            plt.close()
                

    def plotGammasFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.symmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot= []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0])
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                plt.plot(ef_plot, gamma_plot, '-', label = secondParam)

            spin, sublat, orbital, symmetry = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$J_{SC}$ (meV)")
            plt.xlabel(r"$E_{Fermi}$ (meV)")
            plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
            plt.grid()
            plt.xlim(right=-900)
            #plt.xlim(0 , 0.1)
            plt.savefig(f"../Plots/GammaFermi_{spin}_{sublat}_{orbital}_{symmetry}.png")
            plt.close()

    def plotGammasFilling(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.symmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                n_total_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        n_total_plot.append(self.fillingTotal[i]/12.)
                        #gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                plt.plot(n_total_plot, gamma_plot, '-', label = secondParam)

            spin, sublat, orbital, symmetry = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$J_{SC}$ (meV)")
            plt.xlabel(r"$n_{tot}$")
            plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
            plt.grid()
            #plt.xlim(0 , 0.1)
            plt.savefig(f"../Plots/GammaFilling_{spin}_{sublat}_{orbital}_{symmetry}.png")
            plt.close() 

    def plotNnnGammasFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.nnnSymmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot= []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0])
                        gamma_plot.append(np.abs(self.nnnSymmetryGammaDict[key][i]))
                plt.plot(ef_plot, gamma_plot, '-', label = secondParam)

            spin, sublat, orbital, symmetry = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$J_{SC}$ (meV)")
            plt.xlabel(r"$E_{Fermi}$ (meV)")
            plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
            plt.grid()
            #plt.xlim(0 , 0.1)
            plt.savefig(f"../Plots/nnnGammaFermi_{spin}_{sublat}_{orbital}_{symmetry}.png")
            plt.close()

    def plotNnnGammasFilling(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.nnnSymmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                n_total_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        n_total_plot.append(self.fillingTotal[i]/12.)
                        #gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                        gamma_plot.append(np.abs(self.nnnSymmetryGammaDict[key][i]))
                plt.plot(n_total_plot, gamma_plot, '-', label = secondParam)

            spin, sublat, orbital, symmetry = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$J_{SC}$ (meV)")
            plt.xlabel(r"$n_{tot}$")
            plt.ylabel(fr"$\Gamma_{symmetry}$ (meV)")
            plt.grid()
            #plt.xlim(0 , 0.1)
            plt.savefig(f"../Plots/nnnGammaFilling_{spin}_{sublat}_{orbital}_{symmetry}.png")
            plt.close() 

    def plotFillingFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        plt.figure(0)
        for secondParam in secondParamValues:
            ef_plot = []
            n_total_plot = []
            n_chosen_plot_lat1 = []
            n_chosen_plot_lat2 = []

            for i in range(len(self.params)):
                if int(self.params[i][1]) == secondParam:
                    n_total_plot.append(self.fillingTotal[i]/12.)
                    n_chosen_plot_lat1.append(self.filling[(1,1,1)][i]/12.) #key (spin,sublat,orbital)
                    n_chosen_plot_lat2.append(self.filling[(1,2,1)][i]/12.)
                    #gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                    ef_plot.append(self.params[i][0] - self.eMinimal)
            plt.plot(ef_plot, n_chosen_plot_lat1, '-', label = secondParam)
            plt.plot(ef_plot, n_chosen_plot_lat2, '--')
        
        plt.legend(title = r"$U_{Hub}$ (meV)")
        plt.xlabel(r"$E_{Fermi}$ (meV)")
        plt.ylabel(r"$n_{orb}$")
        plt.grid()
        plt.savefig(f"../Plots/FillingFermiOrbital.png")
        plt.close()

        plt.figure(1)
        plt.plot(ef_plot, n_total_plot, '-', label = secondParam)
        plt.legend(title = r"$U_{Hub}$ (meV)")
        plt.xlabel(r"$E_{Fermi}$ (meV)")
        plt.ylabel(r"$n_{tot}$")
        plt.grid()
        plt.savefig(f"../Plots/FillingFermiTotal.png")
        plt.close()


    def plotGammaFermiUnsymmetrized(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        noSymKeys = []
        for spin in (1,2):
            for neighbor in (1,2,3):
                for sublat in (1,2):
                    for orbital in (1,2,3):
                        key = (spin, neighbor, sublat, orbital)
                        noSymKeys.append(key)

        for key in noSymKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot= []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0])
                        gamma_plot.append(np.abs(self.gamma[key][i]))
                plt.plot(ef_plot, gamma_plot, '-', label = secondParam)

            spin, neighbor, sublat, orbital = key
            plt.title(fr's = {-spin + 1.5}, $\alpha$ = {sublat}, l = {orbital}')
            plt.legend(title = r"$J_{SC}$ (meV)")
            plt.xlabel(r"$E_{Fermi}$ (meV)")
            plt.ylabel(fr"$\Gamma_{neighbor}$ (meV)")
            plt.grid()
            plt.xlim(right=-900)
            #plt.xlim(0 , 0.1)
            plt.savefig(f"../Plots/noSymGammaFermi_{spin}_{sublat}_{orbital}_{neighbor}.png")
            plt.close()
       
