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
    eMin = -1053
    #eMin = -1480
    gammaAndFillingPlotter = GammaAndFillingPlotter(runsPath = '/net/ascratch/people/plgjczarnecki/LAO-STO-NNN-corrected', matchPattern = 'RUN_.*', nNeighbors = 3, nNextNeighbors=6, eMinimal=eMin)

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=('e_fermi', 'j_sc_nnn'), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    #gammaAndFillingPlotter.plotGammasFermi()
    #gammaAndFillingPlotter.plotGammasFilling()
    gammaAndFillingPlotter.plotNnnGammasFermi()
    #gammaAndFillingPlotter.plotNnnGammasFilling()
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasTemperature()
    # gammaAndFillingPlotter.plotGammasTemperatureMap()

    #dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadSuperconductingGapMap('/home/jczarnecki/LAO-STO-results/LAO-STO-v0', 'RUN_.*')
    # dispersionPlotter.plotSuperconductingGapMap()
    # Ef_tab = [-1400, -1280]
    # dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadDos(f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{Ef_tab[0]}.0_J_SC_400.0/OutputData/DOS.dat")
    # dos1 = dispersionPlotter.dosDataframe.DOS
    # dispersionPlotter = DispersionPlotter()
    # dispersionPlotter.LoadDos(f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{Ef_tab[1]}.0_J_SC_400.0/OutputData/DOS.dat")
    # dos2 = dispersionPlotter.dosDataframe.DOS
    # E = dispersionPlotter.dosDataframe.E
    # plt.figure(figsize = (6,8), dpi = 400)
    # plt.plot(dos1, E, color="black", linewidth=1, label = "87.5")
    # plt.plot(dos2, E, color="red", linewidth=1, label = "207.5")
    # plt.ylim(bottom=-15, top=15)
    # plt.xlabel(r"DOS")
    # plt.ylabel(r"E (meV)")
    # plt.legend(title = r"$E_{Fermi}$ (meV)")
    # # plt.grid(True)
    # plt.savefig("../Plots/DOS_combined")
    # plt.close()


    #     dispersionPlotter.LoadSuperconductingGap(f'/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{ef}.0_J_SC_400.0/OutputData/SuperconductingGap.dat')
    #     dispersionPlotter.plotSuperconductingGap(postfix = '_J400_'+str(ef), title = fr'$E_{{Fermi}} = {ef - eMin}$ (meV)')

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
    # dispersionPlotter = DispersionPlotter()
    # for ef in [-1280, -1320, -1400]:
    #     dispersionPlotter.LoadSuperconductingGap(
    #         f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{ef}.0_J_SC_200.0/OutputData/SuperconductingGap.dat"
    #     )
    #     dispersionPlotter.plotSuperconductingGap(
    #         postfix=f"_J_100_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
    #     )

    # dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.LoadDos("/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_-1320.0_J_SC_300.0/OutputData/DOS.dat")
    #dispersionPlotter.GetStatistics()
    #dispersionPlotter.shiftEnergies()

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKy", 500, "ky", 0.0, 2, 12
    )

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKx", 500, "kx", 0.0, 2, 12
    )

    # dispersionPlotter.plotFermiCrossection(200, 1.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")

    #dispersionPlotter.plotDos(0.1, "../Plots/DOS.png")
    # for ef in [-1280, -1320, -1400]:
    #     for j_sc in [200, 300]:
    #         dispersionPlotter = DispersionPlotter()
    #         dispersionPlotter.LoadDos(f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0/OutputData/DOS.dat")
    #         dispersionPlotter.plotDos(2, f"../Plots/DOS_Ef{ef}_J{j_sc}.png")
            # dispersionPlotter.LoadSuperconductingGap(
            #     f"/net/ascratch/people/plgjczarnecki/KTO-test/RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0/OutputData/SuperconductingGap.dat"
            # )
            # dispersionPlotter.plotSuperconductingGap(
            #     postfix=f"_J_{j_sc}_Ef_{ef}", title=rf"$E_{{Fermi}} = {ef - eMin}$ (meV)"
            # )


if __name__ == "__main__":
    main()
    # testing()
