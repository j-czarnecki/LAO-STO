import matplotlib.pyplot as plt
import numpy as np
from DataReaderClass import *
from SymmetryResolverClass import *
import seaborn as sns
from scipy.interpolate import griddata
from matplotlib.colors import PowerNorm, Normalize
from matplotlib.cm import ScalarMappable


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
        material: str,
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
        self.material = material
        self.a_tilde = self.__getMaterialsLatticeConstant(self.material)
        self.symmetryKeys: list[tuple[int, int, int, str]] = []
        self.nnnSymmetryKeys: list[tuple[int, int, int, str]] = []
        self.orbitalNameMapping = list[str]
        self.spinSymbolsMapping = list[str]
        self.latticeNameMapping = list[str]
        self.subbandNameMapping = list[str]
        self.maxval = np.float64
        self.efMaxval = np.float64
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
                    self.efMaxval = self.params[i][0]

        if not self.nNextNeighbors == 0:
            for key in self.nnnSymmetryKeys:
                for i in range(len(self.nnnSymmetryGammaDict[key][:])):
                    if np.abs(self.nnnSymmetryGammaDict[key][i]) > self.maxval:
                        self.maxval = np.abs(self.nnnSymmetryGammaDict[key][i])
                        self.efMaxval = self.params[i][0]


        print("Maxval is ", self.maxval)
        print("E_Fermi for maxval is ", self.efMaxval)

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

    def plotGammasFermi(self, eMax: float = None):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))
        gamma_plot = []
        ef_plot = []
        n_total_plot = []

        #cmap = plt.cm.viridis
        #norm = Normalize(vmin=min(secondParamValues), vmax=max(secondParamValues))


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
                        n_total_plot.append(self.__calculateFillingPerCm2(self.fillingTotal[i], self.a_tilde))

                #color = cmap(norm(secondParam))
                #ax1.plot(ef_plot, gamma_plot, label=secondParam, color=color)
                ax1.plot(ef_plot, gamma_plot, label=secondParam)

            band, spin, sublat, orbital, symmetry = key

            ax1.set_ylim(top=1.02 * self.maxval)
            ax1.set_xlim(right = eMax if eMax != None else max(ef_plot))
            ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
            ax1.set_ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )

            ax1.legend(title=r"$J$ (meV)", loc="best")

            # sm = ScalarMappable(cmap=cmap, norm=norm)
            # sm.set_array([])  # Required for ScalarMappable
            # colorbar = fig.colorbar(sm, ax=ax1)
            # colorbar.set_label(r"$J$ (meV)")  # Update label as needed


            # Do this as a last step and trigger plt.draw() so that the ticks are already set in their final form
            plt.draw()
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
            ax2.set_xlabel(r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10)  # Customize units as needed
            plt.savefig(
                f"../Plots/GammaFermi_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotSymmetryRatios(self, eMax: float = None):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        paramsWithoutSym = set(
            ((key[0], key[1], key[2], key[3]) for key in self.symmetryKeys)
        )
        symmetries = set(((key[-1]) for key in self.symmetryKeys))

        self.__setPalette(nColors=len(symmetries))

        colors = ['black', 'red']
        i_color = 0

        ef_plot = []
        n_total_plot = []
        fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
        for key in paramsWithoutSym:
            for secondParam in secondParamValues:
                if key == (1,1,1,1) and secondParam in [50, 100]:
                    #fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
                    gammaSRatio = []
                    gammaPRatio = []
                    gammaDRatio = []
                    ef_plot = []
                    n_total_plot = []
                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            total = 0
                            for sym in symmetries:
                                gammaKey = (*key, sym)
                                total += np.abs(self.symmetryGammaDict[gammaKey][i])
                            gammaKey = (*key, "s")
                            gammaSRatio.append(
                                np.abs(self.symmetryGammaDict[gammaKey][i]) / total * 100
                            )
                            gammaKey = (*key, "p+")
                            gammaPRatio.append(
                                np.abs(self.symmetryGammaDict[gammaKey][i]) / total * 100
                            )
                            gammaKey = (*key, "d+")
                            gammaDRatio.append(
                                np.abs(self.symmetryGammaDict[gammaKey][i]) / total * 100
                            )
                            ef_plot.append(self.params[i][0] - self.eMinimal)
                            n_total_plot.append(self.__calculateFillingPerCm2(self.fillingTotal[i], self.a_tilde))
                    if i_color == 0:
                        ax1.plot(ef_plot, gammaSRatio, label="s", marker=".", markersize=6, color = colors[i_color])
                        ax1.plot(ef_plot, gammaPRatio, label="p+", marker=4, markersize=6, color = colors[i_color])
                        ax1.plot(ef_plot, gammaDRatio, label="d+", marker=5, markersize=6, color = colors[i_color])
                    else:
                        ax1.plot(ef_plot, gammaSRatio, marker=".", markersize=6, color = colors[i_color])
                        ax1.plot(ef_plot, gammaPRatio, marker=4, markersize=6, color = colors[i_color])
                        ax1.plot(ef_plot, gammaDRatio, marker=5, markersize=6, color = colors[i_color])

                    band, spin, sublat, orbital = key
                    # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
                    i_color += 1
                    print(i_color)
        ax1.legend(title=r"$\xi$", loc="center left")
        ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
        ax1.set_ylabel(r"$\Gamma_\xi$ (\%)", labelpad=20)
        ax1.set_ylim(bottom=0, top=100)
        ax1.set_xlim(left = 6, right = eMax if eMax != None else max(ef_plot))

        # Do this as a last step and trigger plt.draw() so that the ticks are already set in their final form
        plt.draw()
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
        ax2.set_xlabel(r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10)  # Customize units as needed
        plt.savefig(
            f"../Plots/AAAGammaSymmetryRatios_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_J_SC_{secondParam}.png"
        )
        plt.close()

    def plotNnnSymmetryRatios(self, eMax: float = None):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        paramsWithoutSym = set(
            ((key[0], key[1], key[2], key[3]) for key in self.nnnSymmetryKeys)
        )
        symmetries = set(((key[-1]) for key in self.nnnSymmetryKeys))

        self.__setPalette(nColors=len(symmetries))

        ef_plot = []
        n_total_plot = []
        for key in paramsWithoutSym:
            for secondParam in secondParamValues:
                fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
                gammaSRatio = []
                gammaPPlusRatio = []
                gammaPMinusRatio = []
                gammaDPlusRatio = []
                gammaDMinusRatio = []
                gammaFRatio = []

                ef_plot = []
                n_total_plot = []
                for i in range(len(self.params)):
                    if int(self.params[i][1]) == secondParam:
                        total = 0
                        for sym in symmetries:
                            gammaKey = (*key, sym)
                            total += np.abs(self.nnnSymmetryGammaDict[gammaKey][i])
                        gammaKey = (*key, "s")
                        gammaSRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "p+")
                        gammaPPlusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "p-")
                        gammaPMinusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "d+")
                        gammaDPlusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "d-")
                        gammaDMinusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "f")
                        gammaFRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        n_total_plot.append(self.__calculateFillingPerCm2(self.fillingTotal[i], self.a_tilde))

                ax1.plot(ef_plot, gammaSRatio, label="s", marker=".", markersize=6)
                ax1.plot(ef_plot, gammaPPlusRatio, label="p+", marker=4, markersize=6)
                ax1.plot(ef_plot, gammaPMinusRatio, label="p-", marker=5, markersize=6)
                ax1.plot(ef_plot, gammaDPlusRatio, label="d+", marker=6, markersize=6)
                ax1.plot(ef_plot, gammaDMinusRatio, label="d-", marker=7, markersize=6)
                ax1.plot(ef_plot, gammaFRatio, label="f", marker="*", markersize=6)

                band, spin, sublat, orbital = key
                # plt.title(fr'$\sigma$ = {-spin + 1.5}, $\alpha$ = {sublat}, l = {self.orbitalNameMapping[orbital - 1]}')
                ax1.legend(title=r"$\xi$", loc="center left")
                ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
                ax1.set_ylabel(r"$\Gamma_\xi$ (\%)", labelpad=20)
                ax1.set_ylim(bottom=0, top=100)
                ax1.set_xlim(left = 6, right = eMax if eMax != None else max(ef_plot))

                # Do this as a last step and trigger plt.draw() so that the ticks are already set in their final form
                plt.draw()
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
                ax2.set_xlabel(r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10)  # Customize units as needed

                plt.savefig(
                    f"../Plots/nnnGammaSymmetryRatios_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_J_SC_{secondParam}.png"
                )
                plt.close()


    ############################################################
    # TODO: For methods below band parameter should be added!!!!!
    ############################################################

    def plotGammasTemperature(self, matchSecondParam: list = None, matchThirdParam: list = None):
        """Here temperature has to be passed as 0-th argument to GammaAndFilling plotter object, j_sc as the 1-st and e_fermi: 2-nd."""
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        cmap = plt.cm.viridis
        norm = Normalize(vmin=min(thirdParamValues) - self.eMinimal, vmax=max(thirdParamValues) - self.eMinimal)
        if matchThirdParam is not None:
            norm = Normalize(vmin=min(matchThirdParam) - self.eMinimal, vmax=max(matchThirdParam) - self.eMinimal)

        for key in self.symmetryKeys:
            fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
            for thirdParam in thirdParamValues:
                if matchThirdParam is not None:
                    if thirdParam not in matchThirdParam:
                        continue
                for secondParam in secondParamValues:
                    if matchSecondParam is not None:
                        if secondParam not in matchSecondParam:
                            continue
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
                    color = cmap(norm(thirdParam - self.eMinimal))
                    ax1.plot(T_plot, gamma_plot, color=color)

            band, spin, sublat, orbital, symmetry = key

            ax1.set_xlabel(r"$T$ (K)")
            ax1.set_ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )

            # Add colorbar
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  # Required for ScalarMappable
            colorbar = fig.colorbar(sm, ax=ax1)
            colorbar.set_label(r"$E_{Fermi}$ (meV)")  # Update label as needed

            plt.savefig(
                f"../Plots/GammaTemperature_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
            )
            plt.close()

    def plotGammasTemperatureMap(self, eMax: float = None):
        """Temperature has to be passed a 0-th argument to this method"""
        # J_SC
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        # E_Fermi
        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        for key in self.symmetryKeys:
            fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)

            for secondParam in secondParamValues:
                T = []
                Ef = []
                Gap = []
                N_total = []
                for thirdParam in thirdParamValues:
                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            T.append(self.params[i][0])
                            Ef.append(self.params[i][2] - self.eMinimal)
                            Gap.append(np.abs(self.symmetryGammaDict[key][i]))
                            N_total.append(self.__calculateFillingPerCm2(self.fillingTotal[i], self.a_tilde))

                T_unique = np.unique(T)  # Find unique T values
                Ef_unique = np.unique(Ef)  # Find unique Ef values
                N_total_unique = np.unique(N_total)
                Ef_grid, T_grid = np.meshgrid(Ef_unique, T_unique)

                points = np.array([Ef, T]).T
                values = Gap
                Gap_grid = griddata(
                    points, values, (Ef_grid, T_grid), method="linear", fill_value=0
                )

                colormesh = ax1.pcolormesh(
                    Ef_grid, T_grid, Gap_grid, cmap="inferno", norm=PowerNorm(gamma=0.8)
                )

                band, spin, sublat, orbital, symmetry = key
                ax1.set_xlabel(rf"$E_{{Fermi}}$ (meV)")
                ax1.set_ylabel(rf"T (K)")
                ax1.set_xlim(left = 50, right = eMax if eMax != None else max(Ef))

                colorbar = fig.colorbar(colormesh, ax=ax1)
                colorbar.set_label(rf"{self.__getNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)")

                plt.draw()
                ax1_ticks = ax1.get_xticks()
                # Plot secondary axis for occupation
                tick_labels = np.interp(
                    ax1_ticks, Ef, N_total
                )  # Interpolate the mapping
                ax2 = ax1.secondary_xaxis("top")
                ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                ax2.set_xticklabels(
                    [f"{val:.2f}" for val in tick_labels]
                )  # Map `n_total_plot` as tick labels
                ax2.set_xlabel(r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10)  # Customize units as needed

                ax1.tick_params(axis='x', direction = 'out')
                ax1.tick_params(axis='y', direction = 'out')
                ax2.tick_params(axis='x', direction = 'out')

                plt.savefig(
                    f"../Plots/GammaTemperatureMap{secondParam}_band{band}_spin{spin}_lat{sublat}_orb{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
                )
                plt.close()

    def plotNnnGammasFermi(self, eMax: float = None):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        gamma_plot = []
        ef_plot = []
        n_total_plot = []
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
                        n_total_plot.append(self.__calculateFillingPerCm2(self.fillingTotal[i], self.a_tilde))

                ax1.plot(ef_plot, gamma_plot, label=secondParam)

            band, spin, sublat, orbital, symmetry = key
            ax1.legend(title=r"$J$ (meV)", loc="best")
            ax1.set_ylim(top=1.02 * self.maxval)
            ax1.set_xlim(left = 0, right = eMax if eMax != None else max(ef_plot))
            ax1.set_xlabel(r"$E_{Fermi}$ (meV)")
            ax1.set_ylabel(
                rf"{self.__getNextNearestNeighborGammaLabel(sublat, orbital, symmetry, spin)} (meV)",
                labelpad=20,
            )
            # Do this as a last step and trigger plt.draw() so that the ticks are already set in their final form
            plt.draw()
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
            ax2.set_xlabel(r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10)  # Customize units as needed

            plt.savefig(
                f"../Plots/nnnGammaFermi_band{band}_spin{spin}_lat{sublat}_{self.orbitalNameMapping[orbital - 1]}_{symmetry}.png"
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

    def __getMaterialsLatticeConstant(self, material: str) -> float:
        """ Calculates lattice constant of hexagonal lattice for a given material."""
        if material == "KTO":
            return np.sqrt(2./3.) * 3.988 # angstroms
        elif material == "STO":
            return np.sqrt(2./3.) * 3.905 # angstroms
        else:
            raise ValueError("Unknown material")

    def __calculateFillingPerCm2(self, filling: float, a_tilde: float) -> float:
        """ Calculates electronic filling in 10^14 cm^-2 """
        return filling / ( 3 * np.sqrt(3.) / 2. * (a_tilde * 1e-8)**2 ) / 1e14

    def __detectOutliers(self, x, y):
        """This method detects too steep gradients for a given curve and classify them as outliers"""
        dy_dx = np.gradient(y, x)
        threshold = 2 * np.std(dy_dx)
        outliers = np.abs(dy_dx) > threshold
        x_out = [x[i] for i in range(len(x)) if outliers[i]]
        y_out = [y[i] for i in range(len(y)) if outliers[i]]
        return x_out, y_out

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """
