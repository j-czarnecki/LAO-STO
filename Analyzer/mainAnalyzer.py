from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *

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


def main():

    #simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_low_U', matchPattern= 'RUN_.*')
    # simulationData.LoadFilling()
    # simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    # simulationData.sortData()

    # symmetryResolver = SymmetryResolver(3)
    # symmetryResolver.CalculateSymmetryGamma(simulationData.gamma)


    gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_low_U', matchPattern= 'RUN_.*', nNeighbors=3, eMinimal = -1053)
    gammaAndFillingPlotter.LoadFilling()
    gammaAndFillingPlotter.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.plotGammasFermi()

    
    # dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos("../OutputData/DOS.dat")
    # dispersionPlotter.GetStatistics()
    # dispersionPlotter.shiftEnergies()


    # dispersionPlotter.plotKx0Crossection('../Plots/DispersionSlice', 100, 2.5, 12)
    # dispersionPlotter.plotKx0Crossection("../Plots/DispersionSliceZoom", 50, 0.5, 12)
    # dispersionPlotter.plotFermiCrossection(10, 0.05, '../Plots/FermiSlice.png')
    # dispersionPlotter.plotDos(100, "../Plots/DOS.png")
    # dispersionPlotter.plotDos(50, "../Plots/DOSzoom.png")

    #plotGammasFermi(simulationData, symmetryResolver)
    #plotGammasFilling(simulationData, symmetryResolver)
    #plotFillingFermi(simulationData, symmetryResolver)

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