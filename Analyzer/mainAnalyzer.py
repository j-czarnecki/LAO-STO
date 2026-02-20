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
from LogParserClass import *
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
    reader = DataReader(runsPath=os.path.join(SCRATCH_PATH, 'KTO-SC', 'KTO-BAND-B_magnitude_B_phi_adaptive_broyden_mix'),
                        matchPattern="RUN_.*",
                        sublattices=3,
                        subbands=1,
                        nBands=18,
                        nAllNeighbours=12,)
    plotter = GammaAndFillingPlotter(material='KTO')
    xKeywords = ('B_magnitude', 'B_phi')
    gammasDf = reader.LoadGamma(xKeywords)
    chargesDf = reader.LoadFilling(xKeywords)

    #gammasDf = gammasDf[gammasDf['B_phi'].isin([0., 30., 60., 90.])]

    plotter.plotGammasTwoParam2d(gammasDf,
                                 chargesDf,
                                 xKeywords,
                                 firstXLabel=r"$|\vec{B}|$ (T)",
                                 plotSecondX=False,
                                 secondXLabel=r"$n$",
                                 legendTitle=r"$\mu$ (meV)",
                                 yUnit=r"($\mu$eV)",
                                 plotRaw=True,
                                 continuousColor=True)

def plotDispersions():
    reader = DataReader(runsPath="/home/czarnecki/LAO-STO/",
                        matchPattern="RUN_.*",
                        sublattices=3,
                        subbands=1,
                        nBands=18,
                        nAllNeighbours=12,)
    dispersionPlotter = DispersionPlotter(plotOutputPath="../Plots")

    #Dispersion plots
    # dispersionDf = reader.LoadDispersion("../OutputData/Energies.dat")
    # dispersionPlotter.GetStatistics(dispersionDf)
    # dispersionPlotter.plotCrossection(dispersionDf, 500, "ky", 0.0, 2, False)
    # dispersionPlotter.plotCrossection(dispersionDf, 500, "kx", 0.0, 2, False)
    # dispersionPlotter.plotFermiCrossection(dispersionDf, 30, 1.5, False)

    # #DOS plots
    # dosDf = reader.LoadDos("../OutputData/DOS.dat")
    # dispersionPlotter.plotDos(dosDf, 400.0, False, 0.1)

    scGapDf = reader.LoadSuperconductingGap("../OutputData/SuperconductingGap.dat")
    dispersionPlotter.plotSuperconductingGap(scGapDf, "", "")
    dispersionPlotter.plotSuperconductingGapAngular(scGapDf, "", "")

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

def analyzeLogs():
    logParser = LogParser()
    #df = logParser.getDivergentChunks("../log.log")
    df = logParser.getDivergentChunks(os.path.join(SCRATCH_PATH,
                                                   "KTO-SC",
                                                   "KTO-B_phi_B_J_SC_chunking_500_eps4",
                                                   "RUN_B_phi_0.0_B_magnitude_5.0_E_Fermi_0.0",
                                                   "log.log"))
    logParser.plotDivergentChunks(df)

def main():
    logger.info("Starting Analyzer")
    #plotGammas()
    plotDispersions()
    #addMissingBandNumber()
    #analyzeLogs()

if __name__ == "__main__":
    main()
