from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *

SCRATCH_PATH = os.getenv('SCRATCH')

def plotGammas():
    eMin = -1053
    # eMin = -1480
    #eMin = 0
    gammaAndFillingPlotter = GammaAndFillingPlotter(
        runsPath=os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC_J_SC_NNN"),
        matchPattern="RUN_.*",
        nNeighbors=3,
        nNextNeighbors=6,
        eMinimal=eMin,
        sublattices=2,
        subbands=1,
        material='STO'
    )

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "j_sc"), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    #gammaAndFillingPlotter.plotGammasFermi()
    gammaAndFillingPlotter.plotSymmetryRatios(eMax = 100)
    #gammaAndFillingPlotter.plotNnnGammasFermi()
    #gammaAndFillingPlotter.plotNnnSymmetryRatios(eMax = 100)
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    #gammaAndFillingPlotter.plotGammasTemperature()
    #gammaAndFillingPlotter.plotGammasTemperatureMap(eMax = 120)


def plotDispersions():
    eMin = -1053
    dispersionPlotter = DispersionPlotter(sublattices=2, subbands=0)
    for j_sc in [50, 75, 100]:
        for ef in [-1026, -966]:
            dispersionPlotter.LoadSuperconductingGap(os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0_J_SC_NNN_75.0", "OutputData", "SuperconductingGap.dat"))
            dispersionPlotter.plotSuperconductingGap(postfix=f"_mixed_Ef{ef}_J_SC_{j_sc}", title=rf"$E_\text{{Fermi}} = {ef - eMin}$~meV")

    # for ef in [-1280, -1320, -1400]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         os.path.join(SCRATCH_PATH, f"RUN_E_Fermi_{ef}.0_J_SC_200.0", "OutputData", "SuperconductingGap.dat")
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # for j_sc in range (370, 380, 10):
    #     for ef in range(50, 130, 10):
    #         try:
    #             dispersionPlotter.LoadDos(
    #                 os.path.join(
    #                     SCRATCH_PATH,
    #                     "KTO-SC",
    #                     "KTO-fit-E_Fermi_100-150meV",
    #                     f"RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0_SUBLATTICES_3_SUBBANDS_2",
    #                     "OutputData",
    #                     "DOS.dat"
    #                 )
    #             )
    #             dispersionPlotter.plotDos(0.8, f"../Plots/DOS_J{j_sc}_Ef_{ef}.png", False, 8e-2)
    #         except:
    #             continue


    # fermiList = [ef for ef in range(50,130, 10)]
    # dirsList = [os.path.join(
    #                     SCRATCH_PATH,
    #                     "KTO-SC",
    #                     "KTO-fit-E_Fermi_100-150meV",
    #                     f"RUN_E_Fermi_{ef}.0_J_SC_370.0_SUBLATTICES_3_SUBBANDS_2",
    #                     "OutputData",
    #                     "DOS.dat"
    #                 ) for ef in fermiList]
    # dispersionPlotter.plotStackedDos(0.8, f"../Plots/DOS_stack.png", True, 8e-2, dirsList, fermiList)



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
    # dispersionPlotter.plotDos(0.7, "../Plots/DOS.png", False, 1e-2)


def main():
    plotGammas()
    plotDispersions()


if __name__ == "__main__":
    main()
