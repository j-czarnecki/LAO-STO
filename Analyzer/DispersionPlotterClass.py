from DataReaderClass import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import seaborn as sns


# TODO: self.lowestEnergy should not be used - all energies should be shown with respect to E_Fermi
# Data reader is in fact not used here, rethink this architecture
class DispersionPlotter(DataReader):

    def __init__(self, sublattices: int):
        DataReader.__init__(self, "./", "xxx", sublattices)
        self.dataLength: int = 0
        self.kPoints1D: int = 0
        self.maxBands: int = 0

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
        self.maxBands = np.max(self.dispersionDataframe.N)
        self.lowestEnergy = np.min(self.dispersionDataframe.E)
        print(f"Data length is {self.dataLength}")
        print(f"Lowest energy is {self.lowestEnergy} (meV)")
        print(f"Number of k-points is {self.kPoints1D}")
        print(f"Number of bands is {self.maxBands}")

    def shiftEnergies(self):
        lowestEnergy = np.min(self.dispersionDataframe.E)
        self.dispersionDataframe.E -= lowestEnergy
        if not self.dosDataframe.empty:
            self.dosDataframe.E -= lowestEnergy

    def plotFirstBrillouinZoneBoundary(self):
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
        plt.plot(
            brillouinZoneVertices[:, 0],
            brillouinZoneVertices[:, 1],
            "--",
            color="black",
            linewidth=0.8,
        )

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

        plotEnergies = np.zeros((self.kPoints1D, nBands), dtype=np.float32)
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
                        self.dispersionDataframe["P_lat1"][i],
                        self.dispersionDataframe["P_lat2"][i],
                        (
                            self.dispersionDataframe["P_lat3"][i]
                            if "P_lat3" in self.dispersionDataframe
                            else 0
                        ),
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
            np.sort(list(set(self.dispersionDataframe[sliceAlong]))),
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
        self.plotFirstBrillouinZoneBoundary()

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
        plt.xlim(-2.5, 2.5)
        plt.ylim(-2.5, 2.5)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.savefig(plotOutputPath)
        plt.close()

    def plotDos(self, eMax: float, plotOutputPath: str):
        plt.figure()
        plt.plot(self.dosDataframe.E, self.dosDataframe.DOS, color="black", linewidth=1)
        plt.xlim(left=-eMax, right=eMax)
        plt.xlabel(r"DOS")
        plt.ylabel(r"E (meV)")
        # plt.grid(True)
        plt.savefig(plotOutputPath)
        plt.close()

    def plotSuperconductingGap(self, postfix: str, title: str):
        plt.figure()
        self.plotFirstBrillouinZoneBoundary()

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

        # print(self.superconductingGapMap[1][1][1])
        for state in range(1, 7):
            plt.figure()
            self.plotFirstBrillouinZoneBoundary()

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
