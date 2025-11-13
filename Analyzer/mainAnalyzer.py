# This file is part of LAO-STO.
#
# Copyright (C) 2025 Julian Czarnecki
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# If you use this code for scientific research, please cite:
# J. Czarnecki et. al.,
# "Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface",
# arXiv:2508.05075 (2025).
# https://arxiv.org/abs/2508.05075

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
        runsPath=os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_tensor_ham_unit_tested"),
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
    gammaAndFillingPlotter.LoadGamma(xKeywords=("e_fermi", "j_sc_tensor"), loadUnfinished=True)
    gammaAndFillingPlotter.sortData()
    gammaAndFillingPlotter.CalculateSymmetryGamma()
    gammaAndFillingPlotter.getMaxvalSymmetrizedGamma()
    gammaAndFillingPlotter.plotGammasTwoParam2d(firstXLabel=r"$\mu$ (meV)",
                                                neighborsToPlot=("nearest", ),
                                                plotSecondX=False,
                                                secondXLabel=r"$n$ (10\textsuperscript{14} cm\textsuperscript{-2})",
                                                legendTitles=(r"$J$ (meV)",),
                                                continuousColor=False,
                                                yUnit=r"($\mu$eV)")
    #gammaAndFillingPlotter.plotFillingFermi()
    # gammaAndFillingPlotter.plotGammasThreeParamCmap(neighborsToPlot=("nearest",),
    #                                                 secondXLabel=r"$n$ (10\textsuperscript{14} cm\textsuperscript{-2})",
    #                                                 colorUnit=r"($\mu$eV)",)


def plotDispersions():
    #eMin = -1053
    dispersionPlotter = DispersionPlotter(sublattices=2, subbands=1)

    dispersionPlotter.LoadDispersion("../OutputData/Energies.dat")
    dispersionPlotter.GetStatistics()
    #dispersionPlotter.shiftEnergies()

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKy", 500, "ky", 0.0, 2, False
    )

    dispersionPlotter.plotCrossection(
        "../Plots/DispersionSliceKx", 500, "kx", 0.0, 2, False
    )

    #dispersionPlotter.plotFermiCrossection(0, 1.5, "../Plots")

    dispersionPlotter.plotFermiCrossection(50, 1.5, "../Plots")

    dispersionPlotter.plotFermiCrossection(100, 1.5, "../Plots")


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

    # dispersionPlotter.plotFermiCrossection(50, 2.0, "../Plots")
    # dispersionPlotter.plotFermiCrossection(150, 2.0, "../Plots/FermiSlice150.png")
    # dispersionPlotter.plotFermiCrossection(200, 2.0, "../Plots/FermiSlice200.png")
    # dispersionPlotter.plotFermiCrossection(500, 1.0, "../Plots/FermiSlice500.png")
    # for ef in (50,100,):
    #     dispersionPlotter.LoadSuperconductingGap(
    #                 os.path.join(
    #                     os.path.join("/home", "czarnecki", "LAO-STO"),
    #                     "OutputData",
    #                     #"Gap_A1",
    #                     f"SuperconductingGap.dat",
    #                 )
    #             )
    #     #dispersionPlotter.plotSuperconductingGapAngular(postfix=f"", title=rf"")
        #dispersionPlotter.plotSuperconductingGap(postfix=f"{ef}", title=rf"")
        # dispersionPlotter.LoadDos(
        #                 os.path.join(
        #                     os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0"),
        #                     "OutputData",
        #                     "DOS.dat",
        #                 )
        #             )
        # dispersionPlotter.plotDos(0.65, f"../Plots/DOS_{ef}.png", False, 1e-2)

        # dispersionPlotter.plotGammaKMap(
        #     inputPath=os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN", f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0"),
        #     postfix=f"{ef}",
        #     neighborsToPlot=("next",),
        #     plotFermiCrossection=True,
        #     eFermi = -60,
        #     dE=2
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
    logger.info("Starting Analyzer")
    plotGammas()
    #plotDispersions()
    #addMissingBandNumber()


if __name__ == "__main__":
    main()
