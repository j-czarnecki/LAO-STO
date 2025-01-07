from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *


def plotGammas():
    eMin = -1053
    # eMin = -1480
    gammaAndFillingPlotter = GammaAndFillingPlotter(
        runsPath="/net/ascratch/people/plgjczarnecki/LAO-STO-NNN-corrected",
        matchPattern="RUN_.*",
        nNeighbors=3,
        nNextNeighbors=6,
        eMinimal=eMin,
        sublattices=3,
    )

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(
        xKeywords=("e_fermi", "j_sc_nnn"), loadUnfinished=True
    )
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    # gammaAndFillingPlotter.plotGammasFermi()
    # gammaAndFillingPlotter.plotGammasFilling()
    gammaAndFillingPlotter.plotNnnGammasFermi()
    gammaAndFillingPlotter.plotNnnGammasFilling()
    gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    gammaAndFillingPlotter.plotFillingFermi()
    gammaAndFillingPlotter.plotGammasJ()
    gammaAndFillingPlotter.plotFillingFermi()
    gammaAndFillingPlotter.plotGammasTemperature()
    gammaAndFillingPlotter.plotGammasTemperatureMap()


def plotDispersions():
    # eMin = -1053
    dispersionPlotter = DispersionPlotter(sublattices=3)
    # for ef in [-1280, -1320, -1400]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{ef}.0_J_SC_200.0/OutputData/SuperconductingGap.dat"
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos(
    #     "/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_-1320.0_J_SC_300.0/OutputData/DOS.dat"
    # )
    dispersionPlotter.GetStatistics()
    dispersionPlotter.shiftEnergies()

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKy", 5000, "ky", 0.0, 2, 18
    )

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKx", 5000, "kx", 0.0, 2, 18
    )

    dispersionPlotter.plotFermiCrossection(100, 1.0, "../Plots/FermiSlice100.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")
    # dispersionPlotter.plotDos(0.1, "../Plots/DOS.png")


def main():
    # plotGammas()
    plotDispersions()


if __name__ == "__main__":
    main()
