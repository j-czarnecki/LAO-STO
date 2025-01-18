from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *

SCRATCH_PATH = os.getenv('SCRATCH')

def plotGammas():
    # eMin = -1053
    # eMin = -1480
    eMin = 0
    gammaAndFillingPlotter = GammaAndFillingPlotter(
        runsPath=os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-3-layers"),
        matchPattern="RUN_.*",
        nNeighbors=3,
        nNextNeighbors=0,
        eMinimal=eMin,
        sublattices=3,
        subbands=0,
    )

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "j_sc"), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    gammaAndFillingPlotter.plotGammasFermi()
    #gammaAndFillingPlotter.plotSymmetryRatios()
    # gammaAndFillingPlotter.plotGammasFilling()
    # gammaAndFillingPlotter.plotNnnGammasFermi()
    # gammaAndFillingPlotter.plotNnnGammasFilling()
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasTemperature()
    # gammaAndFillingPlotter.plotGammasTemperatureMap()


def plotDispersions():
    # eMin = -1053
    dispersionPlotter = DispersionPlotter(sublattices=3, subbands=2)
    # dispersionPlotter.LoadSuperconductingGap("../OutputData/SuperconductingGap.dat")
    # dispersionPlotter.plotSuperconductingGap(postfix="_test", title=rf"test")
    # for ef in [-1280, -1320, -1400]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         os.path.join(SCRATCH_PATH, f"RUN_E_Fermi_{ef}.0_J_SC_200.0", "OutputData", "SuperconductingGap.dat")
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    dispersionPlotter.LoadDos(
        os.path.join(
            SCRATCH_PATH,
            "KTO-SC",
            "KTO-3-layers-2-subbands",
            "RUN_E_Fermi_60.0_J_SC_290.0_SUBLATTICES_3_SUBBANDS_2",
            "OutputData",
            "DOS.dat"
        )
    )
    # dispersionPlotter.GetStatistics()
    # dispersionPlotter.shiftEnergies()

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKy", 200, "ky", 0.0, 2, True
    # )

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKx", 200, "kx", 0.0, 2, True
    # )

    # dispersionPlotter.plotFermiCrossection(100, 1.0, "../Plots/FermiSlice100.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")
    dispersionPlotter.plotDos(1, "../Plots/DOS.png")


def main():
    plotGammas()
    #plotDispersions()


if __name__ == "__main__":
    main()
