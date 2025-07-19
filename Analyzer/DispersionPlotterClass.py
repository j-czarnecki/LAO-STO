from DataReaderClass import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import seaborn as sns
from scipy.signal import convolve
from matplotlib.colors import PowerNorm, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import colorcet as cc
from matplotlib.collections import LineCollection
# With inline inset creation:
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap

# TODO: self.lowestEnergy should not be used - all energies should be shown with respect to E_Fermi
# Data reader is in fact not used here, rethink this architecture
class DispersionPlotter(DataReader):

    def __init__(self, sublattices: int, subbands: int):
        DataReader.__init__(self, "./", "xxx", sublattices, subbands)
        self.dataLength: int = 0
        self.kPoints1D: int = 0
        self.maxBands: int = 0

        self.__initializePlotParams()

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

    def plotFirstBrillouinZoneBoundary(self, ax = None):
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
        if ax == None:
          ax = plt.gca()
        ax.plot(
            brillouinZoneVertices[:, 0],
            brillouinZoneVertices[:, 1],
            "--",
            color="black",
            linewidth=2,
        )

    def plotCrossection(
        self,
        plotOutputPath: str,
        maxEnergy: float,
        sliceAlong: str,
        fixedKVal: float,
        kMax: float,
        isSuperconducting: bool = False,
    ):

        fixedK = ""
        xLabelOnPlot = ""
        if sliceAlong == "kx":
            xLabelOnPlot = r"$k_x~(\tilde{a}^{-1})$"
            fixedK = "ky"
        elif sliceAlong == "ky":
            xLabelOnPlot = r"$k_y~(\tilde{a}^{-1})$"
            fixedK = "kx"

        filteredDispersion = self.dispersionDataframe[
            self.dispersionDataframe[fixedK] == fixedKVal
        ]
        filteredDispersion["zeros"] = np.zeros(len(filteredDispersion))


        minEnergy = filteredDispersion["E"].min() - 0.02 * maxEnergy
        if isSuperconducting:
            minEnergy = -maxEnergy

        lat3 = "P_lat3" if "P_lat3" in self.dispersionDataframe else "zeros"
        colorValuesDict = {"orbital": ["P_yz", "P_zx", "P_xy"],
                     "spin": ["P_up", "zeros", "P_down"],
                     "lattice": ["P_lat1", "P_lat2", lat3],
                     "quasiparticle": ["P_elec", "zeros", "P_hole"],
        }

        groups = filteredDispersion.groupby("N")

        kZoomMax = 0.12
        eZoomMin = 0.1 * minEnergy
        eZoomMax = 3

        for colorKey in colorValuesDict.keys():
            fig = plt.figure(figsize=(7, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.9, top=0.9, bottom=0.2)
            ax = fig.add_subplot(gs[0,0])
            #left, bottom, width, height = [0.7, 0.7, 0.2, 0.2]
            #axin = fig.add_axes([left, bottom, width, height])

            for _, group in groups:
                x = group[sliceAlong].values
                y = group["E"].values
                c = group[colorValuesDict[colorKey]].values
                # Make line segments
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                # For color per segment, average adjacent colors
                segment_colors = 0.5 * (c[:-1] + c[1:])

                lc = LineCollection(segments, colors=segment_colors, linewidths=2)
                ax.add_collection(lc)

                # lc = LineCollection(segments, colors=segment_colors, linewidths=2)
                # axin.add_collection(lc)

            ax.grid(True, linestyle=':')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
            ax.set_xlim(-kMax, kMax)
            ax.set_ylim(bottom=minEnergy, top=maxEnergy)
            ax.set_xlabel(xLabelOnPlot)
            ax.set_ylabel(r"$E$ (meV)")

            # axin.xaxis.set_ticks_position('top')
            # axin.yaxis.set_ticks_position('right')
            # axin.set_xlim(-kZoomMax, kZoomMax)
            # axin.set_ylim(bottom=eZoomMin, top=eZoomMax)
            # axin.set_xticks([-0.1, 0.0, 0.1])
            # axin.set_yticks([0, eZoomMax/2, eZoomMax])
            # axin.patch.set_alpha(0.6)

            plt.savefig(f"{plotOutputPath}_{colorKey}.png")
            plt.close()

    def plotFermiCrossection(self, eFermi: float, dE: float, plotOutputPath: str):

        filteredDispersion = self.dispersionDataframe[
            np.abs(self.dispersionDataframe["E"] - eFermi) < dE
        ]
        filteredDispersion["zeros"] = np.zeros(len(filteredDispersion))

        lat3 = "P_lat3" if "P_lat3" in self.dispersionDataframe else "zeros"
        colorValuesDict = {"orbital": ["P_yz", "P_zx", "P_xy"],
                     "spin": ["P_up", "zeros", "P_down"],
                     "lattice": ["P_lat1", "P_lat2", lat3],
                     "quasiparticle": ["P_elec", "zeros", "P_hole"],
        }
        groups = filteredDispersion.groupby("N")


        for colorKey in colorValuesDict.keys():
            fig = plt.figure(figsize=(5, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.22, right=0.95, top=0.95, bottom=0.2)
            ax = fig.add_subplot(gs[0,0])
            self.plotFirstBrillouinZoneBoundary()
            if colorKey == "orbital":
                left, bottom, width, height = [0.7, 0.72, 0.2, 0.2]
                axin = fig.add_axes([left, bottom, width, height])
                self.__plotRGBLegend(axin)

            for _, group in groups:
                colors = group[colorValuesDict[colorKey]].values
                ax.scatter(group["kx"], group["ky"], marker="o", s=0.6, c=colors)

            ax.grid(True, linestyle=':')
            ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
            ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")
            ax.set_xlim(-2.5, 2.5)
            ax.set_ylim(-2.5, 2.5)
            ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.set_aspect("equal", adjustable="box")
            plt.savefig(os.path.join(plotOutputPath, f"FermiCrossection_{colorKey}_Ef_{eFermi}"))
            plt.close()

    def plotDos(
        self,
        eMax: float,
        plotOutputPath: str,
        addSmearing: bool,
        zeta: float,
        isSingle: bool = True,
        ax=None,
        color="black",
    ):
        if ax is None and isSingle == False:
            ValueError("If isSingle is False, ax must be provided")

        if addSmearing:
            # Define the Lorentzian broadening
            lorentzian = lambda E: (zeta / np.pi) / (E**2 + zeta**2)

            # Create Lorentzian kernel
            energy_range = 2 * eMax
            step = self.dosDataframe.E[1] - self.dosDataframe.E[0]
            kernel_size = len(self.dosDataframe.E)
            kernel = lorentzian(np.linspace(-energy_range, energy_range, kernel_size))
            kernel /= np.trapz(kernel, dx=step)  # Normalize the kernel

            # Perform convolution
            dosSmoothed = convolve(
                self.dosDataframe.DOS, kernel, mode="same", method="fft"
            )

        if isSingle:
            fig = plt.figure(figsize=(7, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.23, right=0.95, top=0.95, bottom=0.23)
            ax = fig.add_subplot(gs[0,0])

        if addSmearing:
            ax.plot(self.dosDataframe.E, dosSmoothed, color=color, linewidth=1.5)
        else:
            ax.plot(
                self.dosDataframe.E, self.dosDataframe.DOS, color=color, linewidth=1.5
            )

        if isSingle:
            plt.xlim(left=-eMax, right=eMax)
            plt.xlabel(r"E (meV)")
            plt.ylabel(r"DOS")
            plt.savefig(plotOutputPath)
            plt.close()

    def plotStackedDos(
        self,
        eMax: float,
        plotOutputPath: str,
        addSmearing: bool,
        zeta: float,
        dosDirsList: list,
        colorParamList: list,
    ):
        """
        Plots subsequent DOSes on top of each other. User should provide a new LoadDos call for each DOS.
        Moreover list of parameters from which color should be deduced has to be provided.
        """
        gs = gridspec.GridSpec(len(dosDirsList), 1, hspace=-0.7, left=0.1, right=0.9, top=0.95, bottom=0.1)  # negative hspace causes overlap
        fig = plt.figure(figsize=(8,10))
        #fig, ax = plt.subplots(figsize=(7, 5), dpi=400, sharex=True)
        cmap = plt.cm.viridis
        norm = Normalize(0, vmax=max(colorParamList))
        axes = []
        for i, dir in enumerate(dosDirsList[::-1]):
            ax = fig.add_subplot(gs[i, 0])
            axes.append(ax)
            self.LoadDos(dir)
            self.dosDataframe["DOS"] = self.dosDataframe["DOS"]
            color = cmap(norm(colorParamList[dosDirsList.index(dir)]))
            self.plotDos(
                eMax,
                plotOutputPath,
                addSmearing,
                zeta,
                isSingle=False,
                ax=ax,
                color=color,
            )
            ax.set_facecolor("none")
            ax.set_yticks([])
            ax.set_xlim(left=-eMax, right=eMax)

            if i != 0:
                ax.spines['top'].set_visible(False)
                ax.tick_params(top=False)
            if i != len(dosDirsList) - 1:
                ax.set_xticks([])
                ax.spines['bottom'].set_visible(False)
            else:
                ax.set_xticks([])
                ax.set_xlabel(r"E (a.u.)")


        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Required for ScalarMappable
        colorbar = fig.colorbar(sm, ax=axes)
        colorbar.set_label(r"$\mu$ (meV)")  # Update label as needed

        fig.text(0.04, 0.5, r"DOS (a.u.)", va='center', rotation='vertical')
        plt.savefig(plotOutputPath)
        plt.close()

    def plotSuperconductingGap(self, postfix: str, title: str):
        fig = plt.figure(figsize=(7, 5), dpi=400)
        # Set up GridSpec (1 row, 1 column, with some spacing)
        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.22, right=0.65, top=1, bottom=0.05)
        ax = fig.add_subplot(gs[0,0])
        self.plotFirstBrillouinZoneBoundary()
        norm = PowerNorm(gamma=1.2, vmin=0, vmax=self.superconductingGapDataframe.gap.max())
        scat = ax.scatter(
            self.superconductingGapDataframe.kx,
            self.superconductingGapDataframe.ky,
            c=self.superconductingGapDataframe.gap,
            s=0.5,
            cmap="inferno",
            norm=norm,
        )
        print("Minimal value of gap is ", self.superconductingGapDataframe.gap.min())

        cax = fig.add_axes([0.67, 0.22, 0.05, 0.6])  # [left, bottom, width 5% of figure width, height 75% of figure height]
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='inferno'), cax=cax, orientation='vertical')
        cbar.set_label(r"$\tilde{\Delta}$ (meV)")


        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))

        ax.set_title(title)
        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")
        ax.set_aspect("equal")

        # plt.grid()
        plt.savefig("../Plots/SuperconductingGap" + postfix + ".png")
        plt.close()

    def plotSuperconductingGapAngular(self, postfix: str = "", title: str = None):
        #self.__setPalette(nColors=self.superconductingGapDataframe['state'].nunique(), palette="tab10")
        stateColors = ["black", "red", "blue", "orange"]
        figAngular, axAngular = plt.subplots(figsize=(7, 5), dpi=400)
        #baxAngular = brokenaxes(ylims=((0,1), (15,50)), hspace=0.1)
        figFourier, axFourier = plt.subplots(figsize=(7, 5), dpi=400)

        for state, group in self.superconductingGapDataframe.groupby("state"):
            dPhi = 0.05                                   # the maximum allowed x–gap

            # --- sort by polar angle ----------------------------------------------------
            angles = np.arctan2(group.ky, group.kx) / np.pi
            order   = np.argsort(angles)               # pandas/Series works fine with np.argsort
            x = angles.iloc[order].to_numpy()          # turn into contiguous NumPy arrays
            y = group.gap.iloc[order].to_numpy()

            # --- split data into “continuous” chunks ------------------------------------
            break_pts = np.where(np.diff(x) > dPhi)[0] + 1   # index AFTER each large jump
            segments  = np.split(np.arange(x.size), break_pts)

            # --- draw each chunk separately so Matplotlib never bridges the gap ---------
            for seg in segments:
                axAngular.plot(
                    x[seg], y[seg],
                    label=f"{int(state)}" if seg is segments[0] else None,   # one legend entry
                    color=stateColors[int(state) - 1]
                )
            # angles = np.arctan2(group.ky, group.kx) / (np.pi)
            # sortedIndices = np.argsort(angles)
            # sortedAngles = angles.iloc[sortedIndices]
            # sortedGaps = group.gap.iloc[sortedIndices]
            # axAngular.plot(
            #     sortedAngles,
            #     sortedGaps,
            #     label=f"{int(state)}",
            # )

            # gapFft = np.fft.fft(sortedGaps) / len(sortedGaps)
            # freqs = np.fft.fftfreq(len(sortedGaps), d=np.mean(np.diff(sortedAngles))) / np.pi
            # axFourier.scatter(freqs, np.abs(gapFft), label=f"{int(state)}")

        axAngular.set_title(title)
        axAngular.legend(title="n", loc="upper right")
        axAngular.set_xlabel(r"$\varphi$~($\pi$)")
        axAngular.set_ylabel(r"$\tilde{\Delta}$~(meV)")
        figAngular.savefig(f"../Plots/SuperconductingGapAngular{postfix}.png")
        plt.close(figAngular)

        axFourier.set_title(title)
        axFourier.legend(title="n", loc="upper right")
        axFourier.set_xlabel(r"f ($\pi^{-1}$)")
        axFourier.set_ylabel(r"$\tilde{\Delta}$~(meV)")
        axFourier.set_xlim(-2, 2)
        figFourier.savefig(f"../Plots/SuperconductingGapFourier{postfix}.png")
        plt.close(figFourier)

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

    def plotGammaKMap(self,
                      inputPath: str,
                      postfix: str = "",
                      neighborsToPlot: tuple[str, ...] = ("nearest", ),
                      plotFermiCrossection: bool = False,
                      eFermi: float = 0.0,
                      dE: float = 0.0):
        rowNames = [r"$\arg(\Gamma)$ ($\pi$) ", r"$|\Gamma|$ (meV)"]
        columnNames = ["yz", "zx", "xy"]
        nRows = len(rowNames)

        columnOffsetsDict = {"nearest": 0, "next": 2}

        cmapPhase = self.__shiftCmap(cc.cm.cyclic_tritanopic_cwrk_40_100_c20, -0.75)
        cmapModule = "Greys"

        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.layerCouplings + 1):
                    for i, neighborName in enumerate(neighborsToPlot):

                        fig, axes = plt.subplots(2,
                                                 3,
                                                 figsize=(12, 9),
                                                 sharex=True,
                                                 sharey=True,
                                                 constrained_layout=False)
                        for ax in axes.flatten():
                            ax.set_xlim(-2.5, 2.5)
                            ax.set_ylim(-2.5, 2.5)

                        for i, name in enumerate(rowNames):
                            axes[i, 0].set_ylabel(r"$k_y~(\tilde{a}^{-1})$")


                        for i, name in enumerate(columnNames):
                            axes[nRows - 1, i].set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                            axes[0, i].set_title(f"{name}")


                        for orbital in range(1, 4):
                            self.LoadGammaMap(
                                os.path.join(
                                    inputPath,
                                    "OutputData",
                                    f"GammaK_orb{orbital}_spin{spin}_layer{sublat}_band{band}.dat",
                                )
                            )
                            kxGrid, kyGrid = np.meshgrid(np.unique(self.gammaKDataFrame.iloc[:, 0]),
                                                        np.unique(self.gammaKDataFrame.iloc[:, 1]))
                            kPoints = np.array([self.gammaKDataFrame.iloc[:, 0], self.gammaKDataFrame.iloc[:, 1]]).T

                            columnRe = 2 + columnOffsetsDict[neighborName]
                            columnIm = 3 + columnOffsetsDict[neighborName]

                            #Phase part
                            ax = axes[0, orbital - 1]
                            grid = griddata(kPoints,
                                                np.arctan2(self.gammaKDataFrame.iloc[:, columnIm], self.gammaKDataFrame.iloc[:, columnRe]) / np.pi,
                                                (kxGrid, kyGrid),
                                                method="linear",
                                                fill_value=0)
                            colormesh = ax.pcolormesh(kxGrid, kyGrid, grid, cmap=cmapPhase, norm=PowerNorm(gamma=1.))
                            ax.set_aspect("equal")
                            self.plotFirstBrillouinZoneBoundary(ax)
                            if plotFermiCrossection and self.dispersionDataframe is not None:
                                filteredDispersion = self.dispersionDataframe[
                                    np.abs(self.dispersionDataframe["E"] - eFermi) < dE
                                ]
                                groups = filteredDispersion.groupby("N")
                                for _, group in groups:
                                    colors = group[["P_yz", "P_zx", "P_xy"]].values
                                    ax.scatter(group["kx"], group["ky"], marker="o", s=0.6, c=colors)
                            if orbital == 3:
                                # Manual colorbar axis (left, bottom, width, height) in figure coords
                                cax = fig.add_axes([1.02, 0.55, 0.015, 0.42])  # adjust as needed
                                cbar = fig.colorbar(colormesh, cax=cax)
                                cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
                                cbar.set_label(rowNames[0])

                            # Module squared
                            ax = axes[1, orbital - 1]

                            grid = griddata(kPoints, np.float64(
                                    np.sqrt(self.gammaKDataFrame.iloc[:, columnRe] ** 2
                                    + self.gammaKDataFrame.iloc[:, columnIm] ** 2)
                                ), (kxGrid, kyGrid), method="linear", fill_value=0)
                            colormesh = ax.pcolormesh(kxGrid, kyGrid, grid, cmap=cmapModule, norm=PowerNorm(gamma=1.5))
                            ax.set_aspect("equal")
                            self.plotFirstBrillouinZoneBoundary(ax)

                            if plotFermiCrossection and self.dispersionDataframe is not None:
                                filteredDispersion = self.dispersionDataframe[
                                    np.abs(self.dispersionDataframe["E"] - eFermi) < dE
                                ]
                                groups = filteredDispersion.groupby("N")
                                for _, group in groups:
                                    colors = group[["P_yz", "P_zx", "P_xy"]].values
                                    ax.scatter(group["kx"], group["ky"], marker="o", s=0.6, c=colors)

                            if orbital == 3:
                                # Manual colorbar axis (left, bottom, width, height) in figure coords
                                cax = fig.add_axes([1.02, 0.03, 0.015, 0.42])  # adjust as needed
                                cbar = fig.colorbar(colormesh, cax=cax)
                                cbar.set_label(rowNames[1])

                        fig.subplots_adjust(wspace=0, hspace=0.1, left=0, right=1, top=1, bottom=0)
                        fig.savefig( f"../Plots/GammaKMap_{neighborName}_spin{spin}_layer{sublat}_band{band}_{postfix}.png",
                                bbox_inches="tight",
                        )
                        plt.close()

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __initializePlotParams(self):
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Computer Modern Roman"
        plt.rcParams["font.sans-serif"] = "Computer Modern Sans serif"
        plt.rcParams["font.monospace"] = "Computer Modern Typewriter"
        plt.rcParams["axes.titlesize"] = 30
        plt.rcParams["axes.labelsize"] = 30
        plt.rcParams["xtick.labelsize"] = 26
        plt.rcParams["ytick.labelsize"] = 26
        plt.rcParams["font.size"] = 26
        plt.rcParams["legend.fontsize"] = 20
        plt.rcParams["legend.title_fontsize"] = 24
        # Optionally, add custom LaTeX preamble
        plt.rcParams["text.latex.preamble"] = (
            r"\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}"
        )

        self.__setPalette()

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

        plt.rcParams["axes.xmargin"] = 0.01

    def __setPalette(self, nColors: int = 3, palette: str = "hsv"):
        # Choose a seaborn palette
        # has to specify number of lines
        self.palette = sns.color_palette(palette, nColors)
        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=self.palette)


    def __plotRGBLegend(self, ax):
        size = 512
        image = np.ones((size, size, 3))

        # Triangle vertices for RGB channels
        v1 = np.array([1.0, 0.0])  # Red corner
        v2 = np.array([0.0, 0.0])  # Green corner
        v3 = np.array([0.5, np.sqrt(3)/2])  # Blue corner (top)

        # Coordinate grid
        x = np.linspace(0, 1, size)
        y = np.linspace(0, np.sqrt(3)/2, size)
        X, Y = np.meshgrid(x, y)

        # For each pixel, compute barycentric coordinates
        def barycentric_coords(x, y):
            # Transformation matrix for barycentric coordinates
            detT = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])
            l1 = ((v2[0] - x) * (v3[1] - y) - (v3[0] - x) * (v2[1] - y)) / detT
            l2 = ((v3[0] - x) * (v1[1] - y) - (v1[0] - x) * (v3[1] - y)) / detT
            l3 = 1.0 - l1 - l2
            return l1, l2, l3

        for i in range(size):
            for j in range(size):
                x_val, y_val = X[i, j], Y[i, j]
                l1, l2, l3 = barycentric_coords(x_val, y_val)
                if (l1 >= 0) and (l2 >= 0) and (l3 >= 0):  # Inside triangle
                    image[i, j, :] = [l1, l2, l3]
        # Create inset
        ax.imshow(image, extent=(0, 1, 0, np.sqrt(3)/2), origin='lower', alpha = 0.7)

        # Triangle border
        triangle = Polygon([v1, v2, v3], closed=True, edgecolor='k', fill=False, lw=0.5)
        ax.add_patch(triangle)

        # Labels
        ax.text(*(v1 + np.array([0.3, 0.5])), r'$yz$', color='black', ha='right', va='top', fontsize=28, rotation=-60)
        ax.text(*(v2 + np.array([-0.3, 0.5])), r'$zx$', color='black', ha='left', va='top', fontsize=28, rotation=60)
        ax.text(*(v3 + np.array([0., -0.1])), r'$xy$', color='black', ha='center', va='bottom', fontsize=28)

        ax.axis('off')

    def __shiftCmap(self, cmap, fraction_shift=0.0, name='shifted'):
        """Shift a colormap cyclically by `fraction_shift` (0.0 to 1.0)."""
        N = 256
        colors = cmap(np.linspace(0, 1, N))
        shift = int(N * fraction_shift)
        colors = np.roll(colors, shift=shift, axis=0)
        return ListedColormap(colors, name=name)
