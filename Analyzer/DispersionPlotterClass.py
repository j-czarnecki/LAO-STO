from DataReaderClass import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import seaborn as sns


# TODO: self.lowestEnergy should not be used - all energies should be shown with respect to E_Fermi
# Data reader is in fact not used here, rethink this architecture
class DispersionPlotter(DataReader):

    def __init__(self):
        DataReader.__init__(self, "./", "xxx")
        self.dataLength = 0
        self.kPoints1D = 0

        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Computer Modern Roman"
        plt.rcParams["font.sans-serif"] = "Computer Modern Sans serif"
        plt.rcParams["font.monospace"] = "Computer Modern Typewriter"
        plt.rcParams["axes.titlesize"] = 24
        plt.rcParams["axes.labelsize"] = 24
        plt.rcParams["xtick.labelsize"] = 20
        plt.rcParams["ytick.labelsize"] = 20
        # Optionally, add custom LaTeX preamble
        plt.rcParams["text.latex.preamble"] = (
            r"\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}"
        )

        # Choose a seaborn palette
        palette = sns.color_palette("hsv", 7)  # has to specify number of lines

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=palette)

        # Set rcParams for tighter layout
        plt.rcParams["figure.autolayout"] = True
        plt.rcParams["figure.constrained_layout.use"] = False
        plt.rcParams["axes.linewidth"] = 1.2

        # Set rcParams to show ticks on both left and right sides
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams["xtick.bottom"] = True
        plt.rcParams["ytick.left"] = True
        plt.rcParams["xtick.top"] = True
        plt.rcParams["ytick.right"] = True

        plt.rcParams["legend.fontsize"] = 12
        plt.rcParams["legend.title_fontsize"] = 14

        plt.rcParams["axes.xmargin"] = 0.01

        print("Initialized DispersionPlotter object")
        print(self.kPoints1D)

    def GetStatistics(self):
        self.dataLength = len(self.dispersionDataframe.N)
        self.kPoints1D = len(set(self.dispersionDataframe.kx))
        self.lowestEnergy = np.min(self.dispersionDataframe.E)

    def shiftEnergies(self):
        lowestEnergy = np.min(self.dispersionDataframe.E)
        print(f"Lowest energy is {lowestEnergy} (meV)")
        self.dispersionDataframe.E -= lowestEnergy
        self.dosDataframe.E -= lowestEnergy

    def plotCrossection(
        self,
        plotOutputPath: str,
        maxEnergy: float,
        sliceAlong: str,
        fixedKVal: float,
        kMax: float,
        nBands: int,
    ):

        fixedK = ""
        xLabelOnPlot = ""
        if sliceAlong == "kx":
            xLabelOnPlot = r"$k_x~(\tilde{a}^{-1})$"
            fixedK = "ky"
        elif sliceAlong == "ky":
            xLabelOnPlot = r"$k_y~(\tilde{a}^{-1})$"
            fixedK = "kx"

        plotEnergies = np.zeros((self.kPoints1D, nBands), dtype=np.float64)
        currentIndex = np.zeros(nBands, dtype=int)
        for i in range(self.dataLength):
            if self.dispersionDataframe[fixedK][i] == fixedKVal:
                bandNo = int(self.dispersionDataframe.N[i]) - 1
                plotEnergies[currentIndex[bandNo], bandNo] = self.dispersionDataframe.E[
                    i
                ]
                currentIndex[bandNo] += 1
                plt.figure(0)
                plt.plot(
                    self.dispersionDataframe[sliceAlong][i],
                    self.dispersionDataframe.E[i],
                    marker=".",
                    markersize=2,
                    color=(
                        self.dispersionDataframe.P_yz[i],
                        self.dispersionDataframe.P_zx[i],
                        self.dispersionDataframe.P_xy[i],
                    ),
                )
                plt.figure(1)
                plt.plot(
                    self.dispersionDataframe[sliceAlong][i],
                    self.dispersionDataframe.E[i],
                    marker=".",
                    markersize=2,
                    color=(
                        self.dispersionDataframe.P_lat1[i],
                        0,
                        self.dispersionDataframe.P_lat2[i],
                    ),
                )
                plt.figure(2)
                plt.plot(
                    self.dispersionDataframe[sliceAlong][i],
                    self.dispersionDataframe.E[i],
                    marker=".",
                    markersize=2,
                    color=(
                        self.dispersionDataframe.P_up[i],
                        0,
                        self.dispersionDataframe.P_down[i],
                    ),
                )

        plt.figure(0)
        plt.xlim(-kMax, kMax)
        # plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy
        plt.ylim(bottom=-0.02 * maxEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        # plt.grid()
        plt.savefig(plotOutputPath + "_orbital.png")
        plt.close()

        plt.figure(1)
        plt.xlim(-kMax, kMax)
        # plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.ylim(bottom=-0.02 * maxEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        # plt.grid()
        plt.savefig(plotOutputPath + "_lattice.png")
        plt.close()

        plt.figure(2)
        plt.xlim(-kMax, kMax)
        # plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.ylim(bottom=-0.02 * maxEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        # plt.grid()
        plt.savefig(plotOutputPath + "_spin.png")
        plt.close()

        plt.figure(3)
        plt.plot(
            np.sort(list(set(self.dispersionDataframe.ky))),
            plotEnergies,
            linewidth=1,
            color="black",
        )
        plt.xlim(-kMax, kMax)
        # plt.ylim(bottom = self.lowestEnergy - 0.02*maxEnergy, top = maxEnergy)
        plt.ylim(bottom=-0.02 * maxEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        # plt.grid()
        plt.savefig(plotOutputPath + "_standard.png")
        plt.close()

    def plotFermiCrossection(self, eFermi: float, dE: float, plotOutputPath: str):
        plt.figure()
        for i in range(len(self.dispersionDataframe.N)):
            if np.abs(self.dispersionDataframe.E[i] - eFermi) < dE:
                plt.plot(
                    self.dispersionDataframe.kx[i],
                    self.dispersionDataframe.ky[i],
                    marker=".",
                    markersize=0.6,
                    color=(
                        self.dispersionDataframe.P_yz[i],
                        self.dispersionDataframe.P_zx[i],
                        self.dispersionDataframe.P_xy[i],
                    ),
                )

        plt.title(r"$E_{Fermi} = $ " + str(eFermi) + " (meV)")
        plt.xlabel(r"$k_x~(\tilde{a}^{-1})$")
        plt.ylabel(r"$k_y~(\tilde{a}^{-1})$")
        plt.xlim(-2.15, 2.15)
        plt.ylim(-2.15, 2.15)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.savefig(plotOutputPath)
        plt.close()

    def plotDos(self, eMax: float, plotOutputPath: str):
        plt.figure()
        plt.plot(self.dosDataframe.DOS, self.dosDataframe.E, color="black", linewidth=1)
        plt.ylim(bottom=0, top=eMax)
        plt.xlabel(r"DOS")
        plt.ylabel(r"E (meV)")
        # plt.grid(True)
        plt.savefig(plotOutputPath)
        plt.close()

    def plotSuperconductingGap(self, postfix: str, title: str):
        brillouinZoneVertices = np.zeros((7, 2))  # One more to close the polygon

        brillouinZoneVertices[:, 0] = np.array(
            [
                4.0 * np.pi / (3 * np.sqrt(3.0)),
                2.0 * np.pi / (3 * np.sqrt(3.0)),
                -2.0 * np.pi / (3 * np.sqrt(3.0)),
                -4.0 * np.pi / (3 * np.sqrt(3.0)),
                -2.0 * np.pi / (3 * np.sqrt(3.0)),
                2.0 * np.pi / (3 * np.sqrt(3.0)),
                4.0 * np.pi / (3 * np.sqrt(3.0)),
            ]
        )

        brillouinZoneVertices[:, 1] = np.array(
            [
                0.0,
                -2.0 * np.pi / 3.0,
                -2.0 * np.pi / 3.0,
                0.0,
                2.0 * np.pi / 3.0,
                2.0 * np.pi / 3.0,
                0.0,
            ]
        )
        plt.figure()
        plt.plot(
            brillouinZoneVertices[:, 0],
            brillouinZoneVertices[:, 1],
            "--",
            color="black",
            linewidth=0.8,
        )

        plt.scatter(
            self.superconductingGapDataframe.kx,
            self.superconductingGapDataframe.ky,
            c=self.superconductingGapDataframe.gap,
            s=0.5,
            cmap="inferno",
        )
        print("Minimal value of gap is ", self.superconductingGapDataframe.gap.min())
        # for i in range(6):
        #     for j in range(len(self.superconductingGapDataframe.kx)):
        #         #if 3 == int(self.superconductingGapDataframe.state[j]):
        #         if True:
        #             plt.scatter(
        #                 self.superconductingGapDataframe.kx[j],
        #                 self.superconductingGapDataframe.ky[j],
        #                 c=self.superconductingGapDataframe.gap[j],
        #                 s=0.5,
        #                 cmap="Reds",
        #             )

        plt.colorbar(label=r"$\tilde{\Delta}$ (meV)")
        ticks = plt.gca().get_xticks()
        plt.xticks(ticks)
        plt.yticks(ticks)
        plt.xlim(-2.5, 2.5)
        plt.ylim(-2.5, 2.5)

        plt.title(title)
        plt.xlabel(r"$k_x~(\tilde{a}^{-1})$")
        plt.ylabel(r"$k_y~(\tilde{a}^{-1})$")
        plt.gca().set_aspect("equal", adjustable="box")
        # plt.grid()
        plt.savefig("../Plots/SuperconductingGap" + postfix + ".png")
        plt.close()

    def plotSuperconductingGapMap(self):
        brillouinZoneVertices = np.zeros((7, 2))  # One more to close the polygon

        brillouinZoneVertices[:, 0] = np.array(
            [
                4.0 * np.pi / (3 * np.sqrt(3.0)),
                2.0 * np.pi / (3 * np.sqrt(3.0)),
                -2.0 * np.pi / (3 * np.sqrt(3.0)),
                -4.0 * np.pi / (3 * np.sqrt(3.0)),
                -2.0 * np.pi / (3 * np.sqrt(3.0)),
                2.0 * np.pi / (3 * np.sqrt(3.0)),
                4.0 * np.pi / (3 * np.sqrt(3.0)),
            ]
        )

        brillouinZoneVertices[:, 1] = np.array(
            [
                0.0,
                -2.0 * np.pi / 3.0,
                -2.0 * np.pi / 3.0,
                0.0,
                2.0 * np.pi / 3.0,
                2.0 * np.pi / 3.0,
                0.0,
            ]
        )

        # print(self.superconductingGapMap[1][1][1])
        for state in range(1, 7):
            plt.figure()
            plt.plot(
                brillouinZoneVertices[:, 0],
                brillouinZoneVertices[:, 1],
                "--",
                color="black",
                linewidth=0.8,
            )
            plt.scatter(
                self.superconductingGapMap[state]["kx"],
                self.superconductingGapMap[state]["ky"],
                c=self.superconductingGapMap[state]["gap"],
                s=0.5,
                cmap="inferno",
            )
            plt.colorbar(label=r"$\tilde{\Delta}$ (meV)")
            ticks = plt.gca().get_xticks()
            plt.xticks(ticks)
            plt.yticks(ticks)
            plt.xlim(-2.5, 2.5)
            plt.ylim(-2.5, 2.5)

            plt.title(rf"n = {state}")
            plt.xlabel(r"$k_x~(\tilde{a}^{-1})$")
            plt.ylabel(r"$k_y~(\tilde{a}^{-1})$")
            plt.gca().set_aspect("equal", adjustable="box")
            # plt.grid()
            plt.savefig(f"../Plots/SuperconductingGapMap_n{state}.png")
            plt.close()
