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
            linewidth=1,
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

        minEnergy = filteredDispersion["E"].min() - 0.02 * maxEnergy
        if isSuperconducting:
            minEnergy = -maxEnergy
        # Projecting orbitals as color
        colors = filteredDispersion[["P_yz", "P_zx", "P_xy"]].values
        plt.figure()
        plt.scatter(
            filteredDispersion[sliceAlong],
            filteredDispersion["E"],
            marker=".",
            s=2,
            c=colors,
        )
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom=minEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        plt.savefig(plotOutputPath + "_orbital.png")
        plt.close()

        # Projecting spin as color
        colors = filteredDispersion[
            ["P_up", "P_up", "P_down"]
        ].values  # middle value for green color is a dummy variable and is set to 0 later
        colors[:, 1] = 0.0
        plt.figure()
        plt.scatter(
            filteredDispersion[sliceAlong],
            filteredDispersion["E"],
            marker=".",
            s=2,
            c=colors,
        )
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom=minEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        plt.savefig(plotOutputPath + "_spin.png")
        plt.close()

        # Projecting lattice as color
        lat3 = "P_lat3" if "P_lat3" in self.dispersionDataframe else "P_lat1"
        colors = filteredDispersion[
            ["P_lat1", "P_lat2", lat3]
        ].values  # if P_lat3 does not exist assign dummy column P_lat1
        if lat3 == "P_lat1":
            colors[:, 2] = 0.0

        plt.figure()
        plt.scatter(
            filteredDispersion[sliceAlong],
            filteredDispersion["E"],
            marker=".",
            s=2,
            c=colors,
        )
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom=minEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        plt.savefig(plotOutputPath + "_lattice.png")
        plt.close()

        # Projecting electron-hole occupation as color
        colors = filteredDispersion[
            ["P_elec", "P_elec", "P_hole"]
        ].values  # middle value for green color is a dummy variable and is set to 0 later
        colors[:, 1] = 0.0
        plt.figure()
        plt.scatter(
            filteredDispersion[sliceAlong],
            filteredDispersion["E"],
            marker=".",
            s=2,
            c=colors,
        )
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom=minEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        plt.savefig(plotOutputPath + "_elec_hole.png")
        plt.close()

        # Plotting all bands with line
        groups = filteredDispersion.groupby("N")
        plt.figure()
        for _, group in groups:
            plt.plot(
                group[sliceAlong],
                group["E"],
                linewidth=1,
                color="black",
            )
        plt.axhline(-50, color="red", linewidth=0.5, linestyle='--')
        plt.axhline(50, color="red", linewidth=0.5, linestyle='--')
        plt.axhline(0, color="blue", linewidth=0.5, linestyle='--')
        plt.axhline(200, color="blue", linewidth=0.5, linestyle='--')
        plt.xlim(-kMax, kMax)
        plt.ylim(bottom=minEnergy, top=maxEnergy)
        plt.xlabel(xLabelOnPlot)
        plt.ylabel(r"$E$ (meV)")
        plt.savefig(plotOutputPath + "_standard.png")
        plt.close()

    def plotFermiCrossection(self, eFermi: float, dE: float, plotOutputPath: str):

        plt.figure()
        self.plotFirstBrillouinZoneBoundary()
        filteredDispersion = self.dispersionDataframe[
            np.abs(self.dispersionDataframe["E"] - eFermi) < dE
        ]
        groups = filteredDispersion.groupby("N")
        for _, group in groups:
            colors = group[["P_yz", "P_zx", "P_xy"]].values
            plt.scatter(group["kx"], group["ky"], marker="o", s=0.6, c=colors)

        plt.xlabel(r"$k_x~(\tilde{a}^{-1})$")
        plt.ylabel(r"$k_y~(\tilde{a}^{-1})$")
        plt.xlim(-2.5, 2.5)
        plt.ylim(-2.5, 2.5)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.savefig(os.path.join(plotOutputPath, f"FermiCrossection_orbital_Ef_{eFermi}"))
        plt.close()

        #Spin projection
        plt.figure()
        self.plotFirstBrillouinZoneBoundary()
        filteredDispersion = self.dispersionDataframe[
            np.abs(self.dispersionDataframe["E"] - eFermi) < dE
        ]
        groups = filteredDispersion.groupby("N")
        for _, group in groups:
            red = group["P_up"].values
            blue = group["P_down"].values
            green = np.zeros_like(red)
            colors = np.stack([red, green, blue], axis=1)
            plt.scatter(group["kx"], group["ky"], marker="o", s=0.6, c=colors)

        plt.xlabel(r"$k_x~(\tilde{a}^{-1})$")
        plt.ylabel(r"$k_y~(\tilde{a}^{-1})$")
        plt.xlim(-2.5, 2.5)
        plt.ylim(-2.5, 2.5)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.savefig(os.path.join(plotOutputPath, f"FermiCrossection_spin_Ef_{eFermi}"))
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
            fig, ax = plt.subplots(figsize=(7, 5), dpi=400)

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
        plt.figure()
        self.plotFirstBrillouinZoneBoundary()

        plt.scatter(
            self.superconductingGapDataframe.kx,
            self.superconductingGapDataframe.ky,
            c=self.superconductingGapDataframe.gap,
            s=0.5,
            cmap="inferno",
            norm=PowerNorm(
                gamma=1.2, vmin=0, vmax=self.superconductingGapDataframe.gap.max()
            ),
        )
        print("Minimal value of gap is ", self.superconductingGapDataframe.gap.min())

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

    def plotGammaKMap(self, inputPath: str, plotNnn: bool = False):
        multiplier = 2.5
        colorbarTitles = [
            r"$Re\left( \Gamma  \right)$ (meV)"
            r"$Im\left( \Gamma  \right)$ (meV)"
            r"$\left| \Gamma  \right|^2$ (meV)"
        ]
        positiveCmap = LinearSegmentedColormap.from_list("white_to_red", ["white", "red"])

        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.layerCouplings + 1):
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


                        #Real part NN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, self.gammaKDataFrame.iloc[:, 2], (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap='bwr', norm=PowerNorm(gamma=1.))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\mathfrak{Re}\left( \Gamma  \right)$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()

                        plt.savefig(
                            f"../Plots/GammaKMapRe_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
                            bbox_inches="tight",
                        )
                        plt.close()


                        # Imaginary part NN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, self.gammaKDataFrame.iloc[:, 3], (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap='bwr', norm=PowerNorm(gamma=1.))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\mathfrak{Im}\left( \Gamma  \right)$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()

                        plt.savefig(
                            f"../Plots/GammaKMapIm_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
                            bbox_inches="tight",
                        )
                        plt.close()


                        # Module squared NN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, np.float64(
                                np.sqrt(self.gammaKDataFrame.iloc[:, 2] ** 2
                                + self.gammaKDataFrame.iloc[:, 3] ** 2)
                            ), (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap=positiveCmap, norm=PowerNorm(gamma=2))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\left| \Gamma  \right|$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()

                        plt.savefig(
                            f"../Plots/GammaKMapModuleSquared_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
                            bbox_inches="tight",
                        )
                        plt.close()

                        #Only plot NNN if needed
                        if not plotNnn:
                            continue

                        # Real part NNN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, self.gammaKDataFrame.iloc[:, 4], (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap='bwr', norm=PowerNorm(gamma=1.))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\mathfrak{Re}\left( \Gamma  \right)$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()
                        plt.savefig(
                            f"../Plots/GammaKMapNnnRe_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
                            bbox_inches="tight",
                        )
                        plt.close()

                        # Imaginary part NNN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, self.gammaKDataFrame.iloc[:, 5], (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap='bwr', norm=PowerNorm(gamma=1.))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\mathfrak{Im}\left( \Gamma  \right)$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()

                        plt.savefig(
                            f"../Plots/GammaKMapNnnIm_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
                            bbox_inches="tight",
                        )
                        plt.close()

                        # Module squared part NNN
                        fig = plt.figure(figsize=(7, 5), dpi=400)
                        # Set up GridSpec (1 row, 1 column, with some spacing)
                        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
                        ax = fig.add_subplot(gs[0,0])

                        nnRealGrid = griddata(kPoints, np.float64(
                                np.sqrt(self.gammaKDataFrame.iloc[:, 4] ** 2
                                + self.gammaKDataFrame.iloc[:, 5] ** 2)
                            ), (kxGrid, kyGrid), method="linear", fill_value=0)
                        colormesh = ax.pcolormesh(kxGrid, kyGrid, nnRealGrid, cmap=positiveCmap, norm=PowerNorm(gamma=2))

                        ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
                        ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

                        ax.set_xlim(-2.5 * multiplier, 2.5 * multiplier)
                        ax.set_ylim(-2.5 * multiplier, 2.5 * multiplier)

                        ax.tick_params(axis="x", direction="out")
                        ax.tick_params(axis="y", direction="out")

                        #ax.set_aspect("equal", adjustable="box")

                        colorbar = fig.colorbar(colormesh, ax=ax)
                        colorbar.set_label(r"$\left| \Gamma  \right|$ (meV)")

                        self.plotFirstBrillouinZoneBoundary()
                        plt.savefig(
                            f"../Plots/GammaKMapNnnModuleSquared_orb{orbital}_spin{spin}_layer{sublat}_band{band}.png",
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

