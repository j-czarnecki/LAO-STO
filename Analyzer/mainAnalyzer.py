from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *


def main():

    # simulationData = DataReader(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_low_U', matchPattern= 'RUN_.*')
    # simulationData.LoadFilling()
    # simulationData.LoadGamma(xKeywords=('e_fermi', 'u_hub'))
    # simulationData.sortData()

    # symmetryResolver = SymmetryResolver(3)
    # symmetryResolver.CalculateSymmetryGamma(simulationData.gamma)

    # gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath= '../../LAO-STO-results/LAO-STO-Hub', matchPattern= 'RUN_.*', nNeighbors=3, nNextNeighbors=0, eMinimal = -1053)
    # gammaAndFillingPlotter = GammaAndFillingPlotter(
    #     runsPath="/net/ascratch/people/plgjczarnecki/KTO-test",
    #     matchPattern="RUN_.*",
    #     nNeighbors=3,
    #     nNextNeighbors=0,
    #     eMinimal=0,
    # )

    # gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    # gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "j_sc"), loadUnfinished=True)
    # gammaAndFillingPlotter.sortData()
    # gammaAndFillingPlotter.CalculateSymmetryGamma()
    # gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    # gammaAndFillingPlotter.plotGammasFermi()
    # gammaAndFillingPlotter.plotGammasFilling()
    # gammaAndFillingPlotter.plotNnnGammasFermi()
    # gammaAndFillingPlotter.plotNnnGammasFilling()
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasTemperature()
    # gammaAndFillingPlotter.plotGammasTemperatureMap()

    # eMin = -1053
    # for ef in [-1000, -975, -950]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         f"/net/ascratch/people/plgjczarnecki/LAO-STO-test/RUN_E_Fermi_{ef}.0_J_SC_NNN_100.0/OutputData/SuperconductingGap.dat"
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    dispersionPlotter = DispersionPlotter()
    dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    dispersionPlotter.LoadDos("../OutputData/DOS.dat")
    dispersionPlotter.GetStatistics()
    dispersionPlotter.shiftEnergies()

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKy", 500, "ky", 0.0, 2, 12
    )

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKx", 500, "kx", 0.0, 2, 12
    )

    # dispersionPlotter.plotFermiCrossection(200, 1.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")

    # dispersionPlotter.plotDos(500, "../Plots/DOS.png")
    # dispersionPlotter.plotDos(30, "../Plots/DOSzoom.png")


if __name__ == "__main__":
    main()
    # testing()
