import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from DataReaderClass import *
from SymmetryResolverClass import *
import seaborn as sns
from scipy.interpolate import griddata
from matplotlib.colors import PowerNorm, Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as ticker
import logging

logger = logging.getLogger(__name__)

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
        self.spinSymbolSingletTripletNameMapping = list[str]
        self.latticeSingletTripletNameMapping = list[str]
        self.maxval = np.float64
        self.efMaxval = np.float64
        self.__initializeSymmetryKeys()
        self.__initializeMapping()
        self.__initializePlotParams()
        logger.info("Initialized GammaAndFillingPlotter")

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


        logger.info(f"Maxval is {self.maxval} at Ef = {self.efMaxval}")

    def plotGammasFermi(self,
                        firstXLabel: str = r"$\mu$ (meV)",
                        plotSecondX: bool = True,
                        secondXLabel: str = r"$n$",
                        neighborsToPlot: tuple[str, ...] = ("nearest",),
                        legendTitles: tuple[str, ...] = (r"$J$ (meV)", ),
                        firstXMax: float = np.inf,
                        firstXShift: float = 0,
                        continuousColor: bool = False):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        gammaYPlot = []
        firstXPlot = []
        secondXPlot = []

        secondXCallback = None
        neighborGammasList = []
        gammaLabelsCallbacks = []
        gammaNeighorhoodLabels = []

        #Assign neighbor gammas
        if "nearest" in neighborsToPlot:
            neighborGammasList.append(self.symmetryGammaDict)
            gammaLabelsCallbacks.append(self.__getNearestNeighborGammaLabel)
            gammaNeighorhoodLabels.append("nearest")
        if "next" in neighborsToPlot:
            neighborGammasList.append(self.nnnSymmetryGammaDict)
            gammaLabelsCallbacks.append(self.__getNextNearestNeighborGammaLabel)
            gammaNeighorhoodLabels.append("next")

        # Pick second axis label
        if plotSecondX:
            if secondXLabel == r"$n$":
                secondXCallback = self.__calculateFillingPerSpinOrbital
            elif secondXLabel == r"$n$ (10\textsuperscript{14} cm\textsuperscript{-2})":
                secondXCallback = self.__calculateFillingPerCm2
            else:
                raise ValueError(f"Unknown secondXLabel: {secondXLabel}")
        else:
            secondXCallback = lambda *args, **kwargs: None

        if continuousColor:
            cmap = plt.cm.cividis
            norm = Normalize(vmin=min(secondParamValues), vmax=max(secondParamValues))

        # Main plotting loop
        for nNeighborhood, gammaDict in enumerate(neighborGammasList):
            for key in self.symmetryKeys:
                fig = plt.figure(figsize=(7, 5), dpi=400)
                # Set up GridSpec (1 row, 1 column, with some spacing)
                if continuousColor:
                    gs = gridspec.GridSpec(1, 1, figure=fig, left=0.25, right=0.95, top=0.75, bottom=0.2)
                else:
                    gs = gridspec.GridSpec(1, 1, figure=fig, left=0.25, right=0.95, top=0.75, bottom=0.2)
                ax1 = fig.add_subplot(gs[0,0])

                for secondParam in secondParamValues:
                    gammaYPlot = []
                    secondXPlot = []
                    firstXPlot = []

                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            gammaYPlot.append(np.abs(gammaDict[key][i]))
                            firstXPlot.append(self.params[i][0] - firstXShift)
                            secondXPlot.append(secondXCallback(self.fillingTotal[i] * 100))

                    if continuousColor:
                        color = cmap(norm(secondParam))
                        ax1.plot(firstXPlot, gammaYPlot, label=secondParam, color=color, linewidth=2)
                    else:
                        ax1.plot(firstXPlot, gammaYPlot, label=secondParam)

                band, spin, sublat, symmetry = key

                ax1.set_ylim(top=1.02 * self.maxval) # Guarantee a single scale for all plots
                ax1.set_xlim(right=firstXMax if firstXMax != np.inf else max(firstXPlot))
                ax1.set_xlabel(firstXLabel)
                ax1.set_ylabel(
                    rf"{gammaLabelsCallbacks[nNeighborhood](sublat, symmetry, spin)} (meV)",
                    labelpad=20,
                )
                ax1.xaxis.set_major_locator(ticker.LinearLocator(5))
                ax1.yaxis.set_major_locator(ticker.LinearLocator(5))
                for mu in (31, 79, 141):
                    ax1.scatter(mu, 0.02, marker='v', s=75, color='deeppink', zorder=10, edgecolors='k', linewidth=1)

                ax1.grid(True, linestyle=':')
                if continuousColor:
                    sm = ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])  # Required for ScalarMappable
                    colorbar = fig.colorbar(sm, ax=ax1)
                    colorbar.set_label(legendTitles[nNeighborhood])  # Update label as needed TODO: this should be variable
                    colorbar.set_ticks([30, 40])
                else:
                    ax1.legend(title=legendTitles[nNeighborhood], loc="best") #TODO: this should be variable

                # Do this as a last step and trigger plt.draw() so that the ticks are already set in their final form
                if plotSecondX:
                    plt.draw()
                    ax1_ticks = ax1.get_xticks()
                    # Plot secondary axis for occupation
                    tick_labels = np.interp(
                        ax1_ticks, firstXPlot, secondXPlot
                    )  # Interpolate the mapping
                    ax2 = ax1.secondary_xaxis("top")
                    ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                    ax2.set_xticklabels(
                        [f"{val:.0f}" for val in tick_labels]
                    )  # Map `n_total_plot` as tick labels
                    ax2.set_xlabel(fr"{secondXLabel} (10 \textsuperscript{{-2}})", labelpad=16)
                plt.savefig(
                    f"../Plots/GammaFermi_{gammaNeighorhoodLabels[nNeighborhood]}_band{band}_spin{spin}_lat{sublat}_{symmetry}.png"
                )
                plt.close()

    def plotGammasTemperature(
        self, matchSecondParam: list = None, matchThirdParam: list = None
    ):
        """Here temperature has to be passed as 0-th argument to GammaAndFilling plotter object, j_sc as the 1-st and e_fermi: 2-nd."""
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        cmap = plt.cm.viridis
        norm = Normalize(
            vmin=min(thirdParamValues) - self.eMinimal,
            vmax=max(thirdParamValues) - self.eMinimal,
        )
        if matchThirdParam is not None:
            norm = Normalize(
                vmin=min(matchThirdParam) - self.eMinimal,
                vmax=max(matchThirdParam) - self.eMinimal,
            )

        for key in self.symmetryKeys:
            fig = plt.figure(figsize=(7, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.15, right=0.95, top=0.95, bottom=0.2)
            ax1 = fig.add_subplot(gs[0,0])
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

            band, spin, sublat, symmetry = key

            ax1.set_xlabel(r"$T$ (K)")
            ax1.set_ylabel(
                rf"{self.__getNearestNeighborGammaLabel(sublat, symmetry, spin)} (meV)",
                labelpad=20,
            )

            # Add colorbar
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  # Required for ScalarMappable
            colorbar = fig.colorbar(sm, ax=ax1)
            colorbar.set_label(r"$\mu$ (meV)")  # Update label as needed

            plt.savefig(
                f"../Plots/GammaTemperature_band{band}_spin{spin}_lat{sublat}_{symmetry}.png"
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
            fig = plt.figure(figsize=(7, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.75, bottom=0.25)
            ax1 = fig.add_subplot(gs[0,0])

            for secondParam in secondParamValues:
                T = []
                Ef = []
                Gap = []
                Ef_T0 = []
                N_total_T0 = []
                for thirdParam in thirdParamValues:
                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            T.append(self.params[i][0])
                            Ef.append(self.params[i][2] - self.eMinimal)
                            Gap.append(np.abs(self.symmetryGammaDict[key][i]))
                            if self.params[i][0] == self.params[0][0]:
                                Ef_T0.append(self.params[i][2] - self.eMinimal)
                                N_total_T0.append(
                                    # self.__calculateFillingPerCm2(
                                    #     self.fillingTotal[i], self.a_tilde
                                    # )
                                    self.fillingTotal[i] / 12 * 100
                                )


                T_unique = np.unique(T)  # Find unique T values
                Ef_unique = np.unique(Ef)  # Find unique Ef values
                Ef_grid, T_grid = np.meshgrid(Ef_unique, T_unique)

                points = np.array([Ef, T]).T
                values = Gap
                Gap_grid = griddata(
                    points, values, (Ef_grid, T_grid), method="linear", fill_value=0
                )

                colormesh = ax1.pcolormesh(
                    Ef_grid, T_grid, Gap_grid, cmap="inferno", norm=PowerNorm(gamma=.8, vmin=0.0, vmax=self.maxval)
                )

                band, spin, sublat, symmetry = key
                ax1.set_xlabel(rf"$\mu$ (meV)")
                ax1.set_ylabel(rf"T (K)")
                ax1.set_xlim(left=0, right=eMax if eMax != None else max(Ef))
                ax1.xaxis.set_major_locator(ticker.MultipleLocator(30))
                ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.3))

                colorbar = fig.colorbar(colormesh, ax=ax1)
                colorbar.set_label(
                    rf"{self.__getNearestNeighborGammaLabel(sublat, symmetry, spin)} (meV)"
                )

                plt.draw()
                ax1_ticks = ax1.get_xticks()
                # Plot secondary axis for occupation
                tick_labels = np.interp(
                    ax1_ticks, Ef_T0, N_total_T0
                )  # Interpolate the mapping
                ax2 = ax1.secondary_xaxis("top")
                ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                ax2.set_xticklabels(
                    [f"{val:.0f}" for val in tick_labels]
                )  # Map `n_total_plot` as tick labels
                # ax2.set_xlabel(
                #     r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10
                # )  # Customize units as needed
                ax2.set_xlabel(r"$n$ (10\textsuperscript{ -2})", labelpad = 16)

                ax1.tick_params(axis="x", direction="out")
                ax1.tick_params(axis="y", direction="out")
                ax2.tick_params(axis="x", direction="out")

                plt.savefig(
                    f"../Plots/GammaTemperatureMap{secondParam}_band{band}_spin{spin}_lat{sublat}_{symmetry}.png"
                )
                plt.close()


    def plotNnnGammasTemperatureMap(self, eMax: float = None):
        """Temperature has to be passed a 0-th argument to this method"""
        # J_SC
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        # E_Fermi
        thirdParamValues = [element[2] for element in self.params]
        thirdParamValues = sorted(list(set(thirdParamValues)))

        for key in self.nnnSymmetryKeys:
            fig = plt.figure(figsize=(7, 5), dpi=400)
            # Set up GridSpec (1 row, 1 column, with some spacing)
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.75, bottom=0.25)
            ax1 = fig.add_subplot(gs[0,0])

            for secondParam in secondParamValues:
                T = []
                Ef = []
                Gap = []
                N_total_T0 = []
                Ef_T0 = []
                for thirdParam in thirdParamValues:
                    for i in range(len(self.params)):
                        if int(self.params[i][1]) == secondParam:
                            T.append(self.params[i][0])
                            Ef.append(self.params[i][2] - self.eMinimal)
                            Gap.append(np.abs(self.nnnSymmetryGammaDict[key][i]))
                            if self.params[i][0] == self.params[0][0]:
                                Ef_T0.append(self.params[i][2] - self.eMinimal)
                                N_total_T0.append(
                                    # self.__calculateFillingPerCm2(
                                    #     self.fillingTotal[i], self.a_tilde
                                    # )
                                    self.fillingTotal[i] / 12 * 100
                                )

                T_unique = np.unique(T)  # Find unique T values
                Ef_unique = np.unique(Ef)  # Find unique Ef values
                Ef_grid, T_grid = np.meshgrid(Ef_unique, T_unique)

                points = np.array([Ef, T]).T
                values = Gap
                Gap_grid = griddata(
                    points, values, (Ef_grid, T_grid), method="linear", fill_value=0
                )

                colormesh = ax1.pcolormesh(
                    Ef_grid, T_grid, Gap_grid, cmap="inferno", norm=PowerNorm(gamma=1.2, vmin=0.0, vmax=self.maxval)
                )

                band, spin, sublat, symmetry = key
                ax1.set_xlabel(rf"$\mu$ (meV)")
                ax1.set_ylabel(rf"T (K)")
                ax1.set_xlim(left=0, right=eMax if eMax != None else max(Ef))
                ax1.xaxis.set_major_locator(ticker.MultipleLocator(30))
                ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))

                colorbar = fig.colorbar(colormesh, ax=ax1)
                colorbar.set_label(
                    rf"{self.__getNextNearestNeighborGammaLabel(sublat, symmetry, spin)} (meV)"
                )

                plt.draw()
                ax1_ticks = ax1.get_xticks()
                # Plot secondary axis for occupation
                tick_labels = np.interp(
                    ax1_ticks, Ef_T0, N_total_T0
                )  # Interpolate the mapping
                ax2 = ax1.secondary_xaxis("top")
                ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                ax2.set_xticklabels(
                    [f"{val:.0f}" for val in tick_labels]
                )  # Map `n_total_plot` as tick labels
                # ax2.set_xlabel(
                #     r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10
                # )  # Customize units as needed
                ax2.set_xlabel(r"$n$ (10 \textsuperscript{-2})", labelpad = 16)

                ax1.tick_params(axis="x", direction="out")
                ax1.tick_params(axis="y", direction="out")
                ax2.tick_params(axis="x", direction="out")
                plt.savefig(
                    f"../Plots/nnnGammaTemperatureMap{secondParam}_band{band}_spin{spin}_lat{sublat}_{symmetry}.png"
                )
                plt.close()

    def plotFillingFermi(self):
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        fig = plt.figure(figsize=(7, 5), dpi=400)
        # Set up GridSpec (1 row, 1 column, with some spacing)
        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.75, bottom=0.2)
        ax1 = fig.add_subplot(gs[0,0])
        ef_plot = []
        n_total_plot = []
        n_total_centimeters = []

        for secondParam in secondParamValues:
            ef_plot = []
            n_total_plot = []

            for i in range(len(self.params)):
                if int(self.params[i][1]) == int(secondParam):
                    n_total_plot.append(self.fillingTotal[i] / 12.0 * 100)
                    ef_plot.append(self.params[i][0] - self.eMinimal)
            ax1.plot(
                ef_plot,
                n_total_plot,
                "-",
                label=format(secondParam, ".1f"),
            )

        ax1.legend(title=r"$U$ (meV)", loc="upper right")
        ax1.set_xlabel(r"$\mu$ (meV)")
        ax1.set_ylabel(r"$n$ (10 \textsuperscript{-2})", labelpad=10)
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax1.grid(True, linestyle=':')
        plt.draw()

        ax1_ticks = ax1.get_yticks()
        tick_labels = []
        for tick in ax1_ticks:
            tick_labels.append(self.__calculateFillingPerCm2(tick / 100 * 12))
        # Plot secondary axis for occupation
        ax2 = ax1.secondary_yaxis("right")
        ax2.set_yticks(ax1_ticks)  # Use the same positions as `ef_plot`
        ax2.set_yticklabels([f"{val:.1f}" for val in tick_labels])  # Map `n_total_plot` as tick labels
        ax2.set_ylabel(r"$n$ (10 \textsuperscript{-14} cm\textsuperscript{-2})", labelpad=45, rotation=270)  # Customize units as needed


        plt.savefig(f"../Plots/FillingFermiTotal.png")
        plt.close()



    def plotSymmetryRatios(self, eMax: float = None):
        """ DEPRECATED """
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        paramsWithoutSym = set(
            ((key[0], key[1], key[2], key[3]) for key in self.symmetryKeys)
        )
        symmetries = set(((key[-1]) for key in self.symmetryKeys))

        self.__setPalette(nColors=len(symmetries) // 2)

        colors = ["black", "red"]
        i_color = 0

        ef_plot = []
        n_total_plot = []
        fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
        for key in paramsWithoutSym:
            for secondParam in secondParamValues:
                # if key == (1, 1, 1, 3) and secondParam in [50, 100]:
                fig, ax1 = plt.subplots(figsize=(7, 5), dpi=400)
                gammaSRatio = []
                gammaD1Ratio = []
                gammaD2Ratio = []
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
                            np.abs(self.symmetryGammaDict[gammaKey][i])
                            / total
                            * 100
                        )
                        gammaKey = (*key, "d_1")
                        gammaD1Ratio.append(
                            np.abs(self.symmetryGammaDict[gammaKey][i])
                            / total
                            * 100
                        )
                        gammaKey = (*key, "d_2")
                        gammaD2Ratio.append(
                            np.abs(self.symmetryGammaDict[gammaKey][i])
                            / total
                            * 100
                        )
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        n_total_plot.append(
                            # self.__calculateFillingPerCm2(
                            #     self.fillingTotal[i], self.a_tilde
                            # )
                            self.fillingTotal[i]/12
                        )

                    band, spin, sublat = key
                    i_color += 1

                ax1.plot(ef_plot, gammaSRatio, label=r"$s$", marker=".", markersize=6)
                ax1.plot(ef_plot, gammaD1Ratio, label=r"$d_1$", marker=6, markersize=6)
                ax1.plot(ef_plot, gammaD2Ratio, label=r"$d_2$", marker=7, markersize=6)

                ax1.legend(title=r"$\xi$", loc="center left")
                ax1.set_xlabel(r"$\mu$ (meV)")
                ax1.set_ylabel(r"$\Gamma^\xi$ (\%)", labelpad=16)
                ax1.set_ylim(bottom=0, top=100)
                ax1.set_xlim(left=6, right=eMax if eMax != None else max(ef_plot))

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
                ax2.set_xlabel(r"$n$", labelpad = 16)
                # ax2.set_xlabel(
                #     r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10
                # )  # Customize units as needed

                plt.savefig(
                    f"../Plots/GammaSymmetryRatios_band{band}_spin{spin}_lat{sublat}_J_SC_{secondParam}.png"
                )
                plt.close()

    def plotNnnSymmetryRatios(self, eMax: float = None):
        """ DEPRECATED """
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))

        paramsWithoutSym = set(
            ((key[0], key[1], key[3]) for key in self.nnnSymmetryKeys)
        )
        symmetries = set(((key[-1]) for key in self.nnnSymmetryKeys))

        self.__setPalette(nColors=len(symmetries) // 2)

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
                        # gammaKey = (*key, "p_1")
                        # gammaPPlusRatio.append(
                        #     np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        # )
                        # gammaKey = (*key, "p_2")
                        # gammaPMinusRatio.append(
                        #     np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        # )
                        gammaKey = (*key, "d_1")
                        gammaDPlusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        gammaKey = (*key, "d_2")
                        gammaDMinusRatio.append(
                            np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        )
                        # gammaKey = (*key, "f")
                        # gammaFRatio.append(
                        #     np.abs(self.nnnSymmetryGammaDict[gammaKey][i]) / total * 100
                        # )
                        ef_plot.append(self.params[i][0] - self.eMinimal)
                        n_total_plot.append(
                            # self.__calculateFillingPerCm2(
                            #     self.fillingTotal[i], self.a_tilde
                            # )
                            self.fillingTotal[i]/12
                        )

                ax1.plot(ef_plot, gammaSRatio, label=r"$s$", marker=".", markersize=6)
                # ax1.plot(ef_plot, gammaPPlusRatio, label="p_1", marker=4, markersize=6)
                # ax1.plot(ef_plot, gammaPMinusRatio, label="p_2", marker=5, markersize=6)
                ax1.plot(ef_plot, gammaDPlusRatio, label=r"$d_1$", marker=6, markersize=6)
                ax1.plot(ef_plot, gammaDMinusRatio, label=r"$d_2$", marker=7, markersize=6)
                # ax1.plot(ef_plot, gammaFRatio, label="f", marker="*", markersize=6)

                band, spin, sublat = key
                ax1.legend(title=r"$\xi$", loc="center left")
                ax1.set_xlabel(r"$\mu$ (meV)")
                ax1.set_ylabel(r"$\Gamma^\xi$ (\%)", labelpad=20)
                ax1.set_ylim(bottom=0, top=100)
                ax1.set_xlim(left=6, right=eMax if eMax != None else max(ef_plot))

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
                # ax2.set_xlabel(
                #     r"$n $ (10\textsuperscript{14} cm\textsuperscript{-2})", labelpad=10
                # )  # Customize units as needed
                ax2.set_xlabel(r"$n$", labelpad = 16)

                plt.savefig(
                    f"../Plots/nnnGammaSymmetryRatios_band{band}_spin{spin}_lat{sublat}_J_SC_{secondParam}.png"
                )
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
            plt.xlabel(r"$\mu$ (meV)")
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
                    for symmetry in self.projector.getSymmetryNames():
                        self.symmetryKeys.append(
                            (band, spin, sublat, symmetry)
                        )

        # Next-to-nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.sublattices + 1):
                    for symmetry in self.projector.getSymmetryNames():
                        self.nnnSymmetryKeys.append(
                            (band, spin, sublat, symmetry)
                        )

    def __initializeMapping(self):
        self.orbitalNameMapping = ["yz", "zx", "xy"]
        self.spinSymbolsMapping = [r"\uparrow", r"\downarrow"]
        self.spinSymbolSingletTripletNameMapping = [
            r"S",
            r"T",
        ]
        self.latticeNameMapping = [rf"Ti_{i}" for i in range(1, self.sublattices + 1)]
        self.subbandNameMapping = [rf"n_{i}" for i in range(1, self.subbands + 1)]

    def __initializePlotParams(self):
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Computer Modern Roman"
        plt.rcParams["font.sans-serif"] = "Computer Modern Sans serif"
        plt.rcParams["font.monospace"] = "Computer Modern Typewriter"
        plt.rcParams["axes.titlesize"] = 36
        plt.rcParams["axes.labelsize"] = 36
        plt.rcParams["xtick.labelsize"] = 30
        plt.rcParams["ytick.labelsize"] = 30
        plt.rcParams["legend.fontsize"] = 20
        plt.rcParams["legend.title_fontsize"] = 24
        # Optionally, add custom LaTeX preamble
        plt.rcParams["text.latex.preamble"] = (
            r"\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}"
        )

        self.__setPalette()

        # Set rcParams for tighter layout
        plt.rcParams["figure.autolayout"] = False
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

    def __setPalette(self, nColors: int = 3, palette: str = "contrast5"):
        # Define available custom palettes
        custom_palettes = {
            "contrast5": ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
            "contrast10": ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd','#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
            "contrastScientific9": ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2','#D55E00', '#CC79A7', '#999999', '#000000', '#FFFFFF'],
            "brb3": ['#000000','#000000', "#E4000F","#E4000F", "#013FF9"],
        }

        if palette in custom_palettes:
            base_palette = custom_palettes[palette]
            # Repeat colors if nColors > palette length
            self.palette = [base_palette[i % len(base_palette)] for i in range(nColors)]
        else:
            # Use seaborn palette
            self.palette = sns.color_palette(palette, nColors)

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=self.palette)

    def __getNearestNeighborGammaLabel(
        self, sublat: int, symmetry: str, spin: int
    ) -> str:
        return rf"$\Gamma_{{\alpha \overline{{\alpha}}}}^{{{symmetry}}}$"

    def __getNextNearestNeighborGammaLabel(
        self, sublat: int, symmetry: str, spin: int
    ) -> str:
        return rf"$\Gamma_{{\alpha \alpha}}^{{{symmetry}}}$"

    def __getMaterialsLatticeConstant(self, material: str) -> float:
        """Calculates lattice constant of hexagonal lattice for a given material."""
        if material == "KTO":
            return np.sqrt(2.0 / 3.0) * 3.988  # angstroms
        elif material == "STO":
            return np.sqrt(2.0 / 3.0) * 3.905  # angstroms
        else:
            raise ValueError("Unknown material")

    def __calculateFillingPerCm2(self, filling: float) -> float:
        """Calculates electronic filling in 10^14 cm^-2"""
        return filling / (3 * np.sqrt(3.0) / 2.0 * (self.a_tilde * 1e-8) ** 2) / 1e14

    def __calculateFillingPerSpinOrbital(self, filling: float) -> float:
        return filling / 12

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
