from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *

SCRATCH_PATH = os.getenv("SCRATCH")

def plotGammas():
    eMin = -1053
    #eMin = -1480
    #eMin = 0
    gammaAndFillingPlotter = GammaAndFillingPlotter(
        runsPath=os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC_J_SC_NNN"),
        # runsPath=os.path.join(
        #     "/home", "jczarnecki", "LAO-STO-results", "LAO-STO-E_Fermi_J_SC_J_SC_NNN"
        # ),
        matchPattern="RUN_.*",
        nNeighbors=3,
        nNextNeighbors=6,
        eMinimal=eMin,
        sublattices=2,
        subbands=1,
        material="STO",
    )

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "j_sc"), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    #gammaAndFillingPlotter.plotGammasFermi(continuousColor=False, eMax=150)
    #gammaAndFillingPlotter.plotGammasSingleTripletFermi(eMax = 150)
    #gammaAndFillingPlotter.plotSymmetryRatios(eMax=100)
    gammaAndFillingPlotter.plotNnnGammasFermi(eMax=150)
    #gammaAndFillingPlotter.plotNnnSymmetryRatios(eMax = 100)
    # gammaAndFillingPlotter.plotNnnGammasSingleTripletFermi(eMax=150)
    # gammaAndFillingPlotter.plotGammaFermiUnsymmetrized()
    #gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasJ()
    # gammaAndFillingPlotter.plotFillingFermi()
    #gammaAndFillingPlotter.plotGammasTemperature()
    #gammaAndFillingPlotter.plotGammasTemperatureMap(eMax = 150)
    #gammaAndFillingPlotter.plotNnnGammasTemperatureMap(eMax = 150)


def plotDispersions():
    # eMin = -1053
    #dispersionPlotter = DispersionPlotter(sublattices=2, subbands=0)
    # for j_sc in [50, 75, 100]:
    #     for ef in [-1026, -966]:
    #         dispersionPlotter.LoadSuperconductingGap(
    #             os.path.join(
    #                 SCRATCH_PATH,
    #                 "STO-SC",
    #                 "LAO-STO-E_Fermi_J_SC_J_SC_NNN",
    #                 f"RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0_J_SC_NNN_75.0",
    #                 "OutputData",
    #                 "SuperconductingGap.dat",
    #             )
    #         )
    #         dispersionPlotter.plotSuperconductingGap(
    #             postfix=f"_mixed_Ef{ef}_J_SC_{j_sc}",
    #             title=rf"$E_\text{{Fermi}} = {ef - eMin}$~meV",
    #         )

    for ef in [-966, -1026]:
        dispersionPlotter.LoadSuperconductingGap(
            os.path.join(SCRATCH_PATH,
                         "STO-SC",
                         "LAO-STO-E_Fermi_J_SC",
                         f"RUN_E_Fermi_{ef}.0_J_SC_110.0",
                         "OutputData",
                         "SuperconductingGap.dat")
        )
        dispersionPlotter.plotSuperconductingGap(
            postfix=f"Ef_{ef}", title=""
        )

    #dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # for j_sc in range(370, 380, 10):
    #     for ef in range(50, 130, 10):
    #         try:
    #             dispersionPlotter.LoadDos(
    #                 os.path.join(
    #                     SCRATCH_PATH,
    #                     "KTO-SC",
    #                     "KTO-fit-E_Fermi_100-150meV",
    #                     f"RUN_E_Fermi_{ef}.0_J_SC_{j_sc}.0_SUBLATTICES_3_SUBBANDS_2",
    #                     "OutputData",
    #                     "DOS.dat",
    #                 )
    #             )
    #             dispersionPlotter.plotDos(
    #                 0.8, f"../Plots/DOS_J{j_sc}_Ef_{ef}.png", False, 8e-2
    #             )
    #         except:
    #             continue

    # fermiList = [ef for ef in range(50,130, 10)]
    # dirsList = [os.path.join(
    #                     SCRATCH_PATH,
    #                     "KTO-SC",
    #                     "KTO-fit-E_Fermi_100-150meV",
    #                     f"RUN_E_Fermi_{ef}.0_J_SC_360.0_SUBLATTICES_3_SUBBANDS_2",
    #                     "OutputData",
    #                     "DOS.dat"
    #                 ) for ef in fermiList]
    # dispersionPlotter.plotStackedDos(0.8, f"../Plots/DOS_stack.png", True, 8e-2, dirsList, fermiList)

    #dispersionPlotter.GetStatistics()
    #dispersionPlotter.shiftEnergies()

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKy", 300, "ky", 0.0, 2, False
    # )

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKx", 300, "kx", 0.0, 2, False
    # )

    #dispersionPlotter.plotFermiCrossection(200, 1.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(100, 2.0, "../Plots/FermiSlice100.png")
    # dispersionPlotter.plotFermiCrossection(150, 2.0, "../Plots/FermiSlice150.png")
    # dispersionPlotter.plotFermiCrossection(200, 2.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")
    dispersionPlotter = DispersionPlotter(sublattices=3, subbands=1)
    dispersionPlotter.LoadSuperconductingGap(
                os.path.join(
                    "/home/jczarnecki/LAO-STO",
                    "OutputData",
                    "SuperconductingGap.dat",
                )
            )
    dispersionPlotter.plotSuperconductingGapAngular(postfix=f"", title=rf"")
    dispersionPlotter.plotSuperconductingGap(postfix=f"", title=rf"")
    dispersionPlotter.LoadDos(
                    os.path.join(
                        "/home",
                        "jczarnecki",
                        "LAO-STO",
                        "OutputData",
                        "DOS.dat",
                    )
                )
    # dispersionPlotter.plotDos(50, "../Plots/DOS.png", False, 1e-2)

    # dispersionPlotter.plotGammaKMap(
    #     inputPath="/home/jczarnecki/LAO-STO-results/LAO-STO-E_Fermi_J_SC/RUN_E_Fermi_-950.0_J_SC_170.0"
    # )


def addMissingBandNumber():
    runsPath=os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC_NNN")
    matchPattern="RUN_.*_J_SC_NNN_100.0"#"RUN_.*"
    directories = [
        dir for dir in os.listdir(runsPath) if re.match(matchPattern, dir)
    ]

    for dir in directories:
        gammaFile = None
        if os.path.exists(os.path.join(runsPath, dir, "OutputData", "Gamma_SC_final.dat")):
            gammaFile = os.path.join(runsPath, dir, "OutputData", "Gamma_SC_final.dat")
        else:
            gammaFile = os.path.join(runsPath, dir, "OutputData", "Gamma_SC_iter.dat")

        chargeFile = None
        if os.path.exists(os.path.join(runsPath, dir, "OutputData", "Charge_dens_final.dat")):
            chargeFile = os.path.join(runsPath, dir, "OutputData", "Charge_dens_final.dat")
        else:
            chargeFile = os.path.join(runsPath, dir, "OutputData", "Charge_dens_iter.dat")

        files = [gammaFile, chargeFile]
        for file in files:
            with open(file, 'r+') as f:
                firstRow = True
                lines = f.readlines()
                f.seek(0)
                for line in lines:
                    if firstRow:
                        newLine = f"#band {line}\n"
                        f.write(newLine)
                        firstRow = False
                    else:
                        line = line.strip()
                        if line:  # skip empty lines
                            new_line = f"1  {line}\n"
                            f.write(new_line)
                        else:
                            f.write("\n")
                f.truncate()

def main():
    plotGammas()
    #plotDispersions()
    #addMissingBandNumber()


if __name__ == "__main__":
    main()
