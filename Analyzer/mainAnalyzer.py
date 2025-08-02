from DataReaderClass import *
from DispersionPlotterClass import *
from SymmetryResolverClass import *
from GammaAndFillingPlotter import *
import logging

SCRATCH_PATH = os.getenv("SCRATCH")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(filename)s:%(funcName)s:%(lineno)s - %(levelname)s: %(message)s",
    handlers=[logging.FileHandler("analyzer.log"),
              logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def plotGammas():
    #eMin = -1053
    #eMin = -1480
    eMin = 0
    gammaAndFillingPlotter = GammaAndFillingPlotter(
        runsPath=os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN_T"),
        # runsPath=os.path.join(
        #     "/home", "jczarnecki", "LAO-STO-results", "LAO-STO-E_Fermi_J_SC_J_SC_NNN"
        # ),
        matchPattern="RUN_.*",
        nNeighbors=3,
        nNextNeighbors=6,
        eMinimal=eMin,
        sublattices=3,
        subbands=1,
        material="KTO",
    )

    gammaAndFillingPlotter.LoadFilling(loadUnfinished=True)
    gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "t", "j_sc_nnn"), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    # gammaAndFillingPlotter.plotGammasFermi(firstXLabel=r"$\mu$ (meV)",
    #                                        neighborsToPlot=("next", ),
    #                                        plotSecondX=False,
    #                                        secondXLabel=r"$n$ (10\textsuperscript{14} cm\textsuperscript{-2})",
    #                                        legendTitles=(r"$J_{nnn}$ (meV)",),
    #                                        continuousColor=False,
    #                                        firstXShift=eMin,
    #                                        firstXMax=10,
    #                                        yMax=100,
    #                                        yUnit=r"($\mu$eV)")
    #gammaAndFillingPlotter.plotFillingFermi()
    gammaAndFillingPlotter.plotGammasThreeParamCmap(neighborsToPlot=("next",),
                                                    secondXLabel=r"$n$ (10\textsuperscript{14} cm\textsuperscript{-2})",)


def plotDispersions():
    eMin = -1053
    dispersionPlotter = DispersionPlotter(sublattices=3, subbands=1)

    dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.GetStatistics()
    #dispersionPlotter.shiftEnergies()

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKy", 300, "ky", 0.0, 2, False
    # )

    # dispersionPlotter.plotCrossection(
    #     "../Plots/DispersionSliceKx", 300, "kx", 0.0, 2, False
    # )

    # efs = [-60, -52, -44, -36]
    # dosDirs = [os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0", "OutputData", "DOS.dat") for ef in efs]

    # dispersionPlotter.plotStackedDos(
    #         eMax=650,
    #         plotOutputPath="../Plots/DOS_stack.png",
    #         addSmearing=False,
    #         zeta=0.0,
    #         dosDirsList=dosDirs,
    #         colorParamList=efs,
    #     )

    # dispersionPlotter.plotFermiCrossection(85, 1.0, "../Plots")
    # dispersionPlotter.plotFermiCrossection(50, 2.0, "../Plots")
    # dispersionPlotter.plotFermiCrossection(150, 2.0, "../Plots/FermiSlice150.png")
    # dispersionPlotter.plotFermiCrossection(200, 2.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")
    for ef in (-60,):
        # dispersionPlotter.LoadSuperconductingGap(
        #             os.path.join(
        #                 os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0"),
        #                 "OutputData",
        #                 f"SuperconductingGap.dat",
        #             )
        #         )
        #dispersionPlotter.plotSuperconductingGapAngular(postfix=f"", title=rf"")
        #dispersionPlotter.plotSuperconductingGap(postfix=f"{ef}", title=rf"")
        # dispersionPlotter.LoadDos(
        #                 os.path.join(
        #                     os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0"),
        #                     "OutputData",
        #                     "DOS.dat",
        #                 )
        #             )
        # dispersionPlotter.plotDos(0.65, f"../Plots/DOS_{ef}.png", False, 1e-2)

        dispersionPlotter.plotGammaKMap(
            inputPath=os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0"),
            postfix=f"{ef}",
            neighborsToPlot=("next",),
            plotFermiCrossection=True,
            eFermi = -60,
            dE=2
        )


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
    logger.info("Starting Analyzer")
    plotGammas()
    #plotDispersions()
    #addMissingBandNumber()


if __name__ == "__main__":
    main()
