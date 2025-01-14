import matplotlib.pyplot as plt
import numpy as np
from DataReaderClass import *
from SymmetryResolverClass import *
import seaborn as sns
from scipy.interpolate import griddata
from matplotlib.colors import PowerNorm


# TODO: this class should be improved to be more general and possibly plot more symmetries of gamma
# self.eMinimal should not be used, as all energies must be calculated with respect to E_Fermi
class GammaAndFillingPlotter(SymmetryResolver):

    def __init__(
        self,
        runsPath: str,
        matchPattern: str,
        nNeighbors: int,
        nNextNeighbors: int,
        eMinimal: float,
        sublattices: int,
        subbands: int,
    ):
        SymmetryResolver.__init__(
            self,
            nNeighbors,
            nNextNeighbors,
            runsPath,
            matchPattern,
            sublattices,
            subbands,
        )
        self.eMinimal = eMinimal
        self.symmetryKeys: list[tuple[int, int, int, str]] = []
        self.nnnSymmetryKeys: list[tuple[int, int, int, str]] = []
        self.orbitalNameMapping = list[str]
        self.spinSymbolsMapping = list[str]
        self.latticeNameMapping = list[str]
        self.subbandNameMapping = list[str]
        self.maxval = np.float64
        self.__initializeSymmetryKeys()
        self.__initializeMapping()
        self.__initializePlotParams()
        print("Initialized GammaAndFillingPlotter")

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def getMaxvalSymmetrizedGamma(self):
        self.maxval = 0.0
        for key in self.symmetryKeys:
            for i in range(len(self.symmetryGammaDict[key][:])):
                if np.abs(self.symmetryGammaDict[key][i]) > self.maxval:
                    self.maxval = np.abs(self.symmetryGammaDict[key][i])

        if not self.nNextNeighbors == 0:
            for key in self.nnnSymmetryKeys:
                for i in range(len(self.nnnSymmetryGammaDict[key][:])):
                    if np.abs(self.nnnSymmetryGammaDict[key][i]) > self.maxval:
                        self.maxval = np.abs(self.nnnSymmetryGammaDict[key][i])

        print("Maxval is ", self.maxval)

    def plotGammasJ(self):
        E_fermi_tab = [-1030, -1040, -1050]
        colors = ["black", "blue", "red"]
        colorIndex = 0
        for key in self.symmetryKeys:
            plt.figure()
            for ef in E_fermi_tab:
                gamma_plot = []
                j_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == ef:
                        j_plot.append(self.params[i][0])
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                plt.plot(
                    j_plot,
                    gamma_plot,
                    "-",
                    label=ef - self.eMinimal,
                    color=colors[colorIndex % 3],
                )
                colorIndex += 1

            band, spin, sublat, orbital, symmetry = key
            # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
            plt.legend(title=r"$E_{Fermi}$ (meV)")
            plt.xlabel(r"$J$ (meV)")
            plt.ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )
            # plt.grid()
            plt.ylim(0, 0.2)
            plt.savefig(
                f"../Plots/GammaJ_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotGammasFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        for key in self.symmetryKeys:
            fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)

            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot = []
                n_total_plot = []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                        n_total_plot.append(self.fillingTotal[i] / 12.0)

                ax1.plot(ef_plot, gamma_plot, label=secondParam)
                ax1_ticks = ax1.get_xticks()

                # Plot secondary axis for occupation
                tick_labels = np.interp(
                    ax1_ticks, ef_plot, n_total_plot
                )  # Interpolate the mapping
                ax2 = ax1.secondary_xaxis("top")
                ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                ax2.set_xticklabels(
                    [f"{val:.2f}" for val in tick_labels]
                )  # Map `n_total_plot` as tick labels
            band, spin, sublat, orbital, symmetry = key

            ax1.legend(title=r"$J$ (meV)", loc="upper left")
            ax1.set_ylim(top=1.02 * self.maxval)
            ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
            ax2.set_xlabel(r"$n_{total}$", labelpad=10)  # Customize units as needed
            ax1.set_ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )
            plt.savefig(
                f"../Plots/GammaFermi_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotGammasFilling(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        for key in self.symmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                n_total_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        n_total_plot.append(self.fillingTotal[i] / 12.0)
                        # gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                        gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                plt.plot(n_total_plot, gamma_plot, "-", label=secondParam)

            band, spin, sublat, orbital, symmetry = key
            # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
            plt.legend(title=r"$J$ (meV)", loc="upper left")
            plt.xlabel(r"$n_{tot}$")
            plt.ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )
            plt.ylim(top=1.02 * self.maxval)
            # plt.grid()
            # plt.xlim(0 , 0.1)
            plt.savefig(
                f"../Plots/GammaFilling_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotSymmetryRatios(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        paramsWithoutSym = set(
            ((key[0], key[1], key[2], key[3]) for key in self.symmetryKeys)
        )
        symmetries = set(((key[-1]) for key in self.symmetryKeys))

        self.__setPalette(nColors=len(symmetries))

        for key in paramsWithoutSym:
            for secondParam in secondParamValues:
                plt.figure(figsize=(7, 5), dpi=400)
                gammaSRatio = []
                gammaPRatio = []
                gammaDRatio = []
                ef_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        total = 0
                        for sym in symmetries:
                            gammaKey = (*key, sym)
                            total += np.abs(self.symmetryGammaDict[gammaKey][i])
                        gammaKey = (*key, "s")
                        gammaSRatio.append(
                            np.abs(self.symmetryGammaDict[gammaKey][i]) / total
                        )
                        gammaKey = (*key, "p+")
                        gammaPRatio.append(
                            np.abs(self.symmetryGammaDict[gammaKey][i]) / total
                        )
                        gammaKey = (*key, "d+")
                        gammaDRatio.append(
                            np.abs(self.symmetryGammaDict[gammaKey][i]) / total
                        )
                        ef_plot.append(self.params[i][0] - self.eMinimal)

                plt.plot(ef_plot, gammaSRatio, label="s", marker=".", markersize=10)
                plt.plot(ef_plot, gammaPRatio, label="p+", marker=4, markersize=10)
                plt.plot(ef_plot, gammaDRatio, label="d+", marker=5, markersize=10)

                band, spin, sublat, orbital = key
                # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
                plt.legend(title=r"$\xi$", loc="center left")
                plt.xlabel(r"$E_{Fermi}$ (meV)")
                plt.ylabel(r"$\Gamma_\xi$", labelpad=20)
                plt.ylim(bottom=0, top=1)
                plt.savefig(
                    f"../Plots/GammaSymmetryRatios_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_J_SC_{secondParam}.png"
                )
                plt.close()

    ############################################################
    # TODO: For methods below band parameter should be added!!!!!
    ############################################################

    def plotGammasTemperature(self):
        """Here temperature has to be passed as 0-th argument to GammaAndFilling plotter object, j_sc as the 1-st and e_fermi: 2-nd."""
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        thirdParamIterator = 0
        linestyleTab = ["-", "--", ":"]

        for key in self.symmetryKeys:
            plt.figure()
            for thirdParam in thirdParamValues:
                for secondParam in secondParamValues:
                    gamma_plot = []
                    T_plot = []

                    for i in range(len(self.params)):
                        # TODO: change those int(), they may lead to trouble in future
                        if (
                            int(self.params[i][1]) == secondParam
                            and int(self.params[i][2]) == thirdParam
                        ):
                            T_plot.append(self.params[i][0])
                            gamma_plot.append(np.abs(self.symmetryGammaDict[key][i]))
                    plt.plot(T_plot, gamma_plot, linestyle="-", label=secondParam)
                thirdParamIterator += 1
                band, spin, sublat, orbital, symmetry = key
                # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
                plt.title(rf"$E_{{Fermi}} = {thirdParam - self.eMinimal}$ meV", pad=10)
                plt.legend(title=r"$J$ (meV)")
                plt.xlabel(r"$T$ (K)")
                plt.ylabel(
                    rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                    labelpad=20,
                )
                plt.ylim(bottom=-0.02 * self.maxval, top=1.02 * self.maxval)
                plt.xlim(0, 0.6)
                # plt.grid()
                # plt.xlim(0 , 0.1)
                plt.savefig(
                    f"../Plots/GammaTemperatureEf{thirdParam}_{spin}_{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
                )
                plt.close()

    def plotGammasTemperatureMap(self):
        """Temperature has to be passed a 0-th argument to this method"""
        # J_SC
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        # E_Fermi
        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        for key in self.symmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                T = []
                Ef = []
                Gap = []
                for thirdParam in thirdParamValues:
                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            T.append(self.params[i][0])
                            Ef.append(self.params[i][2] - self.eMinimal)
                            Gap.append(np.abs(self.symmetryGammaDict[key][i]))

                T_unique = np.unique(T)  # Find unique T values
                Ef_unique = np.unique(Ef)  # Find unique Ef values
                Ef_grid, T_grid = np.meshgrid(Ef_unique, T_unique)

                points = np.array([Ef, T]).T
                values = Gap
                Gap_grid = griddata(
                    points, values, (Ef_grid, T_grid), method="linear", fill_value=0
                )

                # Gap_grid = np.zeros_like(T_grid)

                # for i in range(len(T)):
                #     row = np.where(T_unique == T[i])[0][0]
                #     col = np.where(Ef_unique == Ef[i])[0][0]
                #     Gap_grid[row, col] = Gap[i]

                plt.pcolormesh(
                    Ef_grid, T_grid, Gap_grid, cmap="inferno", norm=PowerNorm(gamma=0.6)
                )

                spin, sublat, orbital, symmetry = key
                plt.title(rf"$J = {secondParam}$ meV", pad=10)
                plt.xlabel(rf"$E_{{Fermi}}$ (meV)")
                plt.ylabel(rf"T (K)")
                plt.colorbar(
                    label=rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)"
                )
                # plt.ylabel(fr"$\Gamma_{{{self.latticeNameMapping[sublat - 1]} {self.latticeNameMapping[sublat % 2]},{self.orbitalNameMapping[orbital - 1]},{symmetry}}}^{{{self.spinSymbolsMapping[spin - 1]} {self.spinSymbolsMapping[spin % 2]}}}$ (meV)",
                #           labelpad= 20)
                # plt.ylim(bottom = -0.02*self.maxval, top = 1.02*self.maxval)
                # plt.xlim(0, 0.6)
                # plt.grid()
                # plt.xlim(0 , 0.1)
                plt.savefig(
                    f"../Plots/GammaTemperatureMap{secondParam}_{spin}_{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
                )
                plt.close()

    def plotNnnGammasFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.nnnSymmetryKeys:
            fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot = []
                n_total_plot = []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        gamma_plot.append(np.abs(self.nnnSymmetryGammaDict[key][i]))
                        n_total_plot.append(self.fillingTotal[i] / 12.0)

                ax1.plot(ef_plot, gamma_plot, label=secondParam)
                ax1_ticks = ax1.get_xticks()

                # Plot secondary axis for occupation
                tick_labels = np.interp(
                    ax1_ticks, ef_plot, n_total_plot
                )  # Interpolate the mapping
                ax2 = ax1.secondary_xaxis("top")
                ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                ax2.set_xticklabels(
                    [f"{val:.2f}" for val in tick_labels]
                )  # Map `n_total_plot` as tick labels

            band, spin, sublat, orbital, symmetry = key
            ax1.legend(title=r"$J_{nnn}$ (meV)", loc="upper left")
            ax1.set_ylim(top=1.02 * self.maxval)
            ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
            ax2.set_xlabel(r"$n_{total}$", labelpad=10)  # Customize units as needed
            ax1.set_ylabel(
                rf"{self.__getNextNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )
            plt.savefig(
                f"../Plots/nnnGammaFermi_{spin}_{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotNnnGammasFilling(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        for key in self.nnnSymmetryKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                n_total_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        n_total_plot.append(self.fillingTotal[i] / 12.0)
                        # gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                        gamma_plot.append(np.abs(self.nnnSymmetryGammaDict[key][i]))
                plt.plot(n_total_plot, gamma_plot, "-", label=secondParam)

            band, spin, sublat, orbital, symmetry = key
            # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
            plt.legend(title=r"$J_{nnn}$ (meV)", loc="upper right")
            plt.xlabel(r"$n_{tot}$")
            plt.ylabel(
                rf"{self.__getNextNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)"
            )
            plt.ylim(top=1.02 * self.maxval)
            # plt.grid()
            # plt.xlim(0 , 0.1)
            plt.savefig(
                f"../Plots/nnnGammaFilling_{spin}_{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotFillingFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        colors = ["black", "blue", "red", "magenta"]
        colorIndex = 0
        plt.figure(0)
        for secondParam in secondParamValues:
            print(secondParam)
            ef_plot = []
            n_total_plot = []
            n_chosen_plot_lat1 = []
            n_chosen_plot_lat2 = []

            for i in range(len(self.params)):
                if int(self.params[i][1]) == int(secondParam):
                    n_total_plot.append(self.fillingTotal[i] / 12.0)
                    n_chosen_plot_lat1.append(
                        self.filling[(1, 1, 1)][i] / 12.0
                    )  # key (spin,sublat,orbital)
                    n_chosen_plot_lat2.append(self.filling[(1, 2, 1)][i] / 12.0)
                    # gamma_plot.append(np.abs(symmetryResolver.symmetryGammaDict[key][i]))
                    ef_plot.append(self.params[i][0] - self.eMinimal)
            plt.plot(
                ef_plot,
                n_chosen_plot_lat1,
                "-",
                label=format(secondParam, ".1f"),
                color=self.palette[colorIndex],
            )
            plt.plot(ef_plot, n_chosen_plot_lat2, "--", color=self.palette[colorIndex])
            colorIndex += 1

        plt.legend(title=r"$U$ (meV)")
        plt.xlabel(r"$E_{Fermi}$ (meV)")
        plt.ylabel(r"$n_{orb}$")
        # plt.grid()
        plt.savefig(f"../Plots/FillingFermiOrbital.png")
        plt.close()

        plt.figure(1)
        plt.plot(ef_plot, n_total_plot, "-", label=secondParam)
        plt.legend(title=r"$U$ (meV)")
        plt.xlabel(r"$E_{Fermi}$ (meV)")
        plt.ylabel(r"$n_{tot}$")
        plt.grid()
        plt.savefig(f"../Plots/FillingFermiTotal.png")
        plt.close()

    def plotGammaFermiUnsymmetrized(self):
        """DEPRECATED"""
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        noSymKeys = []
        for spin in (1, 2):
            for neighbor in (1, 2, 3):
                for sublat in (1, 2):
                    for orbital in (1, 2, 3):
                        key = (spin, neighbor, sublat, orbital)
                        noSymKeys.append(key)

        for key in noSymKeys:
            plt.figure()
            for secondParam in secondParamValues:
                gamma_plot = []
                ef_plot = []

                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        gamma_plot.append(np.abs(self.gamma[key][i]))
                plt.plot(ef_plot, gamma_plot, "-", label=secondParam)

            spin, neighbor, sublat, orbital = key
            plt.title(
                rf"$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}"
            )
            plt.legend(title=r"$J_{SC}$ (meV)")
            plt.xlabel(r"$E_{Fermi}$ (meV)")
            plt.ylabel(rf"$\Gamma_{neighbor}$ (meV)")
            # plt.grid()
            plt.xlim(right=-900)
            # plt.xlim(0 , 0.1)
            plt.savefig(
                f"../Plots/noSymGammaFermi_{spin}_{sublat}_{self.orbitalNameMapping[orbital - 1]}_{neighbor}.png"
            )
            plt.close()

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __initializeSymmetryKeys(self):
        # Nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.layerCouplings + 1):
                    for orbital in range(1, 4):
                        for symmetry in ("s", "p+", "d+"):
                            self.symmetryKeys.append(
                                (band, spin, sublat, orbital, symmetry)
                            )
        # Next-to-nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.sublattices + 1):
                    for orbital in range(1, 4):
                        for symmetry in ("s", "p+", "d+", "f", "d-", "p-"):
                            self.nnnSymmetryKeys.append(
                                (band, spin, sublat, orbital, symmetry)
                            )

    def __initializeMapping(self):
        self.orbitalNameMapping = ["yz", "zx", "xy"]
        self.spinSymbolsMapping = [r"\uparrow", r"\downarrow"]
        self.latticeNameMapping = [rf"Ti_{i}" for i in range(1, self.sublattices + 1)]
        self.subbandNameMapping = [rf"n_{i}" for i in range(1, self.subbands + 1)]

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

    def __getNearestNeighborGammaLabel(
        self, sublat: int, orbital: int, symmetry: str, spin: int
    ) -> str:
        latSubstring = f"{self.latticeNameMapping[sublat // 2]} {self.latticeNameMapping[(sublat // 2) + sublat % 2 - (sublat + 1) % 2]}"
        spinSubstring = (
            f"{self.spinSymbolsMapping[spin - 1]} {self.spinSymbolsMapping[spin % 2]}"
        )
        return rf"$\Gamma_{{{latSubstring},{self.orbitalNameMapping[orbital - 1]},{symmetry}}}^{{{spinSubstring}}}$"

    def __getNextNearestNeighborGammaLabel(
        self, sublat: int, orbital: int, symmetry: str, spin: int
    ) -> str:
        latSubstring = f"{self.latticeNameMapping[sublat - 1]} {self.latticeNameMapping[sublat - 1]}"
        spinSubstring = (
            f"{self.spinSymbolsMapping[spin - 1]} {self.spinSymbolsMapping[spin % 2]}"
        )
        return rf"$\Gamma_{{{latSubstring},{self.orbitalNameMapping[orbital - 1]},{symmetry}}}^{{{spinSubstring}}}$"

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """
