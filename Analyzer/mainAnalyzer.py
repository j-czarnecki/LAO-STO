from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *

def testing():
    symmetryResolver = SymmetryResolver(3)

    testGammaDict = {}
    for spin in (1,2):
        for sublat in (1,2):
            for orbital in (1,2,3):
                for neighbor in (1,2,3):
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


    gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath= '../../LAO-STO-results/LAO-STO', matchPattern= 'RUN_.*', nNeighbors=3, nNextNeighbors=6, eMinimal = 0)
    gammaAndFillingPlotter.LoadFilling()
    gammaAndFillingPlotter.LoadGamma(xKeywords=('e_fermi', 'j_sc'))
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.plotGammasFermi()
    gammaAndFillingPlotter.plotGammasFilling()
    gammaAndFillingPlotter.plotNnnGammasFermi()
    gammaAndFillingPlotter.plotNnnGammasFilling()
    gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    #gammaAndFillingPlotter.plotGammasJ()
    #gammaAndFillingPlotter.plotFillingFermi()


    # dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos("../OutputData/DOS.dat")
    # dispersionPlotter.GetStatistics()
    # #dispersionPlotter.shiftEnergies()


    # dispersionPlotter.plotCrossection('../Plots/DispersionSliceKy', 200, 'ky', 0.,  2, 12)
    # dispersionPlotter.plotCrossection("../Plots/DispersionSliceZoomKy", 100, 'ky', 0., 0.5, 12)
    
    # dispersionPlotter.plotCrossection('../Plots/DispersionSliceKx', 200, 'kx', 0.,  2, 12)
    # dispersionPlotter.plotCrossection("../Plots/DispersionSliceZoomKx", 100, 'kx', 0., 0.5, 12)
    
    # #dispersionPlotter.plotFermiCrossection(10, 0.05, '../Plots/FermiSlice10.png')
    # #dispersionPlotter.plotFermiCrossection(40, 0.05, '../Plots/FermiSlice40.png')
    # dispersionPlotter.plotFermiCrossection(0, 0.05, '../Plots/FermiSlice0.png')
    # dispersionPlotter.plotFermiCrossection(100, 1, '../Plots/FermiSlice100.png')
    # dispersionPlotter.plotFermiCrossection(150, 0.05, '../Plots/FermiSlice150.png')
    # #dispersionPlotter.plotFermiCrossection(160, 0.05, '../Plots/FermiSlice160.png')
    # dispersionPlotter.plotFermiCrossection(200, 0.05, '../Plots/FermiSlice200.png')

    # dispersionPlotter.plotDos(200, "../Plots/DOS.png")
    # #dispersionPlotter.plotDos(10, "../Plots/DOSzoom.png")



if __name__ == "__main__":
    main()
    #testing()