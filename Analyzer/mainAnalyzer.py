from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *


def testing():
    symmetryResolver = SymmetryResolver(3)

    testGammaDict = {}
    for spin in (1, 2):
        for sublat in (1, 2):
            for orbital in (1, 2, 3):
                for neighbor in (1, 2, 3):
                    key = (spin, neighbor, sublat, orbital)
                    testGammaDict[key] = [x * neighbor for x in range(1, 3)]

    print(testGammaDict)
    symmetryResolver.CalculateSymmetryGamma(testGammaDict)
    print(symmetryResolver.symmetryGammaDict[(1, 1, 1, "s")])
    print(symmetryResolver.symmetryGammaDict[(1, 1, 1, "p")])
    print(symmetryResolver.symmetryGammaDict[(1, 1, 1, "d")])


def main():

    # simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_low_U', matchPattern= 'RUN_.*')
    # simulationData.LoadFilling()
    # simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    # simulationData.sortData()

    # symmetryResolver = SymmetryResolver(3)
    # symmetryResolver.CalculateSymmetryGamma(simulationData.gamma)

    # gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath= '../../LAO-STO-results/LAO-STO-Hub', matchPattern= 'RUN_.*', nNeighbors=3, nNextNeighbors=0, eMinimal = -1053)
    gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath = '/net/ascratch/people/plgjczarnecki/KTO-test', matchPattern = 'RUN_.*', nNeighbors = 3, nNextNeighbors=0, eMinimal=0)

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=('e_fermi', 'j_sc'), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    gammaAndFillingPlotter.plotGammasFermi()
    gammaAndFillingPlotter.plotGammasFilling()
    # gammaAndFillingPlotter.plotNnnGammasFermi()
    # gammaAndFillingPlotter.plotNnnGammasFilling()
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasTemperature()
    # gammaAndFillingPlotter.plotGammasTemperatureMap()

    #dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadSuperconductingGapMap('/home/jczarnecki/LAO-STO-results/LAO-STO-v0', 'RUN_.*')
    # dispersionPlotter.plotSuperconductingGapMap()

    # Ef_tab = [-1000.0, -995.0, -990.0, -985.0, -980.0, -975.0, -970.0, -965.0, -960.0, -955.0, -950.0]

    # for ef in Ef_tab:
    #     dispersionPlotter = DispersionPlotter()
    #     dispersionPlotter.LoadSuperconductingGap(f'/home/jczarnecki/LAO-STO-results/LAO-STO-v0/RUN_E_Fermi_{ef}/OutputData/SuperconductingGap.dat')
    #     dispersionPlotter.plotSuperconductingGap(postfix = '_J170_'+str(ef), title = fr'$E_{{Fermi}} = {ef}$ (meV)')

    # dispersionPlotter.LoadSuperconductingGap(
    #     "/home/jczarnecki/Downloads/LAO-STO-Romberg-test/RUN_E_Fermi_-950.0_J_SC_NNN_100.0/OutputData/SuperconductingGap.dat"
    # )
    # dispersionPlotter.LoadSuperconductingGap(
    #     "/home/jczarnecki/LAO-STO/OutputData/SuperconductingGap.dat"
    # )
    # dispersionPlotter.plotSuperconductingGap(
    #     postfix="_test", title=r"$E_{Fermi} = 53$ (meV)"
    # )

    # eMin = -1053
    # for ef in [-1000, -975, -950]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         f"/net/ascratch/people/plgjczarnecki/LAO-STO-test/RUN_E_Fermi_{ef}.0_J_SC_NNN_100.0/OutputData/SuperconductingGap.dat"
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos("../OutputData/DOS.dat")
    # dispersionPlotter.GetStatistics()
    # dispersionPlotter.shiftEnergies()

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKy", 150, "ky", 0.0, 2, 12
    # )
    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceZoomKy", 30, "ky", 0.0, 0.3, 12
    # )

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKx", 150, "kx", 0.0, 2, 12
    # )
    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceZoomKx", 20, "kx", 0.0, 0.1, 12
    # )

    # dispersionPlotter.plotFermiCrossection(10, 0.05, "../Plots/FermiSlice10.png")
    # dispersionPlotter.plotFermiCrossection(20, 0.05, "../Plots/FermiSlice20.png")
    # dispersionPlotter.plotFermiCrossection(30, 0.05, "../Plots/FermiSlice30.png")
    # dispersionPlotter.plotFermiCrossection(60, 0.05, "../Plots/FermiSlice60.png")

    # dispersionPlotter.plotFermiCrossection(0, 0.5, "../Plots/FermiSlice0.png")
    # dispersionPlotter.plotFermiCrossection(100, 0.5, "../Plots/FermiSlice100.png")
    # dispersionPlotter.plotFermiCrossection(150, 0.5, "../Plots/FermiSlice150.png")
    # dispersionPlotter.plotFermiCrossection(160, 0.5, "../Plots/FermiSlice160.png")
    # dispersionPlotter.plotFermiCrossection(200, 0.5, "../Plots/FermiSlice200.png")

    # dispersionPlotter.plotDos(150, "../Plots/DOS.png")
    # dispersionPlotter.plotDos(30, "../Plots/DOSzoom.png")


if __name__ == "__main__":
    main()
    # testing()
