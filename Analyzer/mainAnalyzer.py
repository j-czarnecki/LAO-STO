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


    # gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath= '../../LAO-STO-results/LAO-STO-Hub', matchPattern= 'RUN_.*', nNeighbors=3, nNextNeighbors=0, eMinimal = -1053)
    
    # gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    # gammaAndFillingPlotter.LoadGamma(xKeywords=('e_fermi', 'j_sc'), loadUnfinished=True)
    # gammaAndFillingPlotter.sortData()
    # gammaAndFillingPlotter.CalculateSymmetryGamma()
    # gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    # gammaAndFillingPlotter.plotGammasFermi()
    #gammaAndFillingPlotter.plotGammasFilling()
    # gammaAndFillingPlotter.plotNnnGammasFermi()
    # gammaAndFillingPlotter.plotNnnGammasFilling()
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasTemperature()


    dispersionPlotter = DispersionPlotter()
    dispersionPlotter.LoadSuperconductingGap('../OutputData/SuperconductingGap_NN_J175.dat')
    dispersionPlotter.plotSuperconductingGap(postfix = 'NN_EF1000_J175', title = r'$J = 175$ (meV)')
    
    dispersionPlotter.LoadSuperconductingGap('../OutputData/SuperconductingGap_NN_J150.dat')
    dispersionPlotter.plotSuperconductingGap(postfix = 'NN_EF1000_J160', title = r'$J = 160$ (meV)')
    
    dispersionPlotter.LoadSuperconductingGap('../OutputData/SuperconductingGap_NNN_J175.dat')
    dispersionPlotter.plotSuperconductingGap(postfix = 'NNN_EF1000_J175', title = r'$J_{nnn} = 175$ (meV)')

    dispersionPlotter.LoadSuperconductingGap('../OutputData/SuperconductingGap_NNN_J150.dat')
    dispersionPlotter.plotSuperconductingGap(postfix = 'NNN_EF1000_J125', title = r'$J_{nnn} = 125$ (meV)')

    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos("../OutputData/DOS.dat")
    # dispersionPlotter.GetStatistics()
    # dispersionPlotter.shiftEnergies()


    # dispersionPlotter.plotCrossection('../Plots/DispersionSliceKy', 150, 'ky', 0.,  2, 12)
    # dispersionPlotter.plotCrossection("../Plots/DispersionSliceZoomKy", 30, 'ky', 0., 0.3, 12)
    
    # dispersionPlotter.plotCrossection('../Plots/DispersionSliceKx', 150, 'kx', 0.,  2, 12)
    # dispersionPlotter.plotCrossection("../Plots/DispersionSliceZoomKx", 30, 'kx', 0., 0.3, 12)

    #dispersionPlotter.plotFermiCrossection(3, 0.05, '../Plots/FermiSlice3.png')
    #dispersionPlotter.plotFermiCrossection(13, 0.05, '../Plots/FermiSlice13.png')
    #dispersionPlotter.plotFermiCrossection(23, 0.05, '../Plots/FermiSlice23.png')


    # dispersionPlotter.plotFermiCrossection(0, 0.05, '../Plots/FermiSlice0.png')
    # dispersionPlotter.plotFermiCrossection(100, 1, '../Plots/FermiSlice100.png')
    # dispersionPlotter.plotFermiCrossection(150, 0.05, '../Plots/FermiSlice150.png')
    # #dispersionPlotter.plotFermiCrossection(160, 0.05, '../Plots/FermiSlice160.png')
    # dispersionPlotter.plotFermiCrossection(200, 0.05, '../Plots/FermiSlice200.png')

    #dispersionPlotter.plotDos(150, "../Plots/DOS.png")
    #dispersionPlotter.plotDos(30, "../Plots/DOSzoom.png")



if __name__ == "__main__":
    main()
    #testing()