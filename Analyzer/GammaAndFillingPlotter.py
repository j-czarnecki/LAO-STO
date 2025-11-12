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
        self.symmetryKeys: dict[str, list[tuple[int, int, int, int, str]]]= {"nearest": [], "next": []}
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
        for key in self.symmetryKeys["nearest"]:
            for i in range(len(self.symmetryGammaDict[key][:])):
                if np.abs(self.symmetryGammaDict[key][i]) > self.maxval:
                    self.maxval = np.abs(self.symmetryGammaDict[key][i])
                    self.efMaxval = self.params[i][0]

        if not self.nNextNeighbors == 0:
            for key in self.symmetryKeys["next"]:
                for i in range(len(self.nnnSymmetryGammaDict[key][:])):
                    if np.abs(self.nnnSymmetryGammaDict[key][i]) > self.maxval:
                        self.maxval = np.abs(self.nnnSymmetryGammaDict[key][i])
                        self.efMaxval = self.params[i][0]


        logger.info(f"Maxval is {self.maxval} at Ef = {self.efMaxval}")

    def plotGammasTwoParam2d(self,
                             firstXLabel: str = r"$\mu$ (meV)",
                             plotSecondX: bool = True,
                             secondXLabel: str = r"$n$",
                             neighborsToPlot: tuple[str, ...] = ("nearest",),
                             legendTitles: tuple[str, ...] = (r"$J$ (meV)", ),
                             firstXMax: float = np.inf,
                             firstXShift: float = 0,
                             yMax: float = np.inf,
                             yUnit: str = "(meV)",
                             continuousColor: bool = False):
        """
        Plots a 2D curve of symmetrized Gammas as a function of argument X,
        where X is the first parameter specified in self.LoadGammas(xKeywords=(X, Y)).
        The second parameter, Y, specifies how many curves, corresponding to different values of Y, are plotted.

        Parameters:
        firstXLabel: str
            Label of the X axis.
        plotSecondX: bool
            If True, a second X axis is plotted.
            It corresponds to carrier concentration and is only sensible for chemical potential/Fermi energy
            being the first argument.
        secondXLabel: str
            Label of the second X axis.
        neighborsToPlot: tuple[str, ...]
            List of neighbors to plot.
            Supported: "nearest", "next".
        legendTitles: tuple[str, ...]
            List of legend titles for nearest/next neighbors. Order must be kept the same as in neighborsToPlot.
            It corresponds to the name of Y parameter.
        firstXMax: float
            Maximum value of the first X axis.
        firstXShift: float
            Shift of the first X axis.
        yMax: float
            Maximum value of the Y axis.
            This is not Y parameter passed to xKeywords.
        yUnit: str
            Unit of the Y axis.
            Supported: '(eV)', '(meV)', '($\\mu$eV) <- Without one backslash'
            This is not Y parameter passed to xKeywords.
        continuousColor: bool
            If True, consecutive curves with different Y values will have a gradually changing color.
            Colorbar will be also plotted.
        """
        secondParamValues = [element[1] for element in self.params]
        secondParamValues = sorted(list(set(secondParamValues)))
        self.__setPalette(nColors=len(set(secondParamValues)))

        gammaYPlot = []
        firstXPlot = []
        secondXPlot = []

        secondXCallback = None
        neighborGammasList = []
        neighborKeys = []
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

        # Set energy units and multipliers
        yMultiplier = 1
        if yUnit == r"(meV)":
            yMultiplier = 1
        elif yUnit == r"($\mu$eV)":
            yMultiplier = 1e3
        elif yUnit == r"(eV)":
            yMultiplier = 1e-3
        else:
            raise ValueError(f"Unknown yUnit: {yUnit}")

        if continuousColor:
            cmap = plt.cm.cividis
            norm = Normalize(vmin=min(secondParamValues), vmax=max(secondParamValues))

        # Main plotting loop
        for nNeighborhood, gammaDict in enumerate(neighborGammasList):
            for key in self.symmetryKeys[gammaNeighorhoodLabels[nNeighborhood]]:
                fig = plt.figure(figsize=(7, 5), dpi=100)
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
                            gammaYPlot.append(np.abs(gammaDict[key][i]) * yMultiplier)
                            firstXPlot.append(self.params[i][0] - firstXShift)
                            secondXPlot.append(secondXCallback(self.fillingTotal[i]))
                            #secondXPlot.append(secondXCallback(self.fillingTotal[i] * 100))

                    if continuousColor:
                        color = cmap(norm(secondParam))
                        ax1.plot(firstXPlot, gammaYPlot, label=secondParam, color=color, linewidth=2)
                    else:
                        ax1.plot(firstXPlot, gammaYPlot, label=secondParam)

                band, spin1, spin2, sublat, symmetry = key

                ax1.set_ylim(bottom=0, top=1.02 * self.maxval * yMultiplier if yMax == np.inf else yMax) # Guarantee a single scale for all plots
                #ax1.set_xlim(right=firstXMax if firstXMax != np.inf else max(firstXPlot))
                ax1.set_xlabel(firstXLabel)
                ax1.set_ylabel(
                    rf"{gammaLabelsCallbacks[nNeighborhood](sublat, symmetry, spin1, spin2)}" + yUnit,
                    labelpad=20,
                )
                # ax1.xaxis.set_major_locator(ticker.LinearLocator(5))
                ax1.yaxis.set_major_locator(ticker.LinearLocator(4))
                ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))

                #for mu in (31, 79, 141):
                    #ax1.scatter(mu, 0.02, marker='v', s=75, color='deeppink', zorder=10, edgecolors='k', linewidth=1)

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
                        [f"{val:.1f}" for val in tick_labels]
                    )  # Map `n_total_plot` as tick labels
                    ax2.set_xlabel(fr"{secondXLabel}", labelpad=16)
                    #ax2.set_xlabel(fr"{secondXLabel} (10 \textsuperscript{{-2}})", labelpad=16)
                plt.savefig(
                    f"../Plots/Gamma2d_{gammaNeighorhoodLabels[nNeighborhood]}_band{band}_spin1{spin1}_spin2{spin2}_lat{sublat}_{symmetry}.png"
                )
                plt.close()

    def plotGammasThreeParamCmap(self,
                                 firstXLabel: str = r"$\mu$ (meV)",
                                 plotSecondX: bool = True,
                                 secondXLabel: str = r"$n$",
                                 neighborsToPlot: tuple[str, ...] = ("nearest",),
                                 firstXMax: float = np.inf,
                                 firstXShift: float = 0,
                                 yMax: float = np.inf,
                                 yUnit: str = "(K)",
                                 colorMax: float = np.inf,
                                 colorUnit: str = "(meV)"):
        """
        Plots a 2D colormap of symmetrized Gammas as a function of arguments X and Y,
        where X is the first parameter specified in self.LoadGammas(xKeywords=(X, Y, Z))
        and Y is the second one.
        The third parameter, Z, specifies how many maps, corresponding to different values of Z, are plotted.

        Parameters:
        firstXLabel: str
            Label of the X axis.
        plotSecondX: bool
            If True, a second X axis is plotted.
            It corresponds to carrier concentration and is only sensible for chemical potential/Fermi energy
            being the first argument.
        secondXLabel: str
            Label of the second X axis.
        neighborsToPlot: tuple[str, ...]
            List of neighbors to plot.
            Supported: "nearest", "next".
        firstXMax: float
            Maximum value of the first X axis.
        firstXShift: float
            Shift of the first X axis.
        yMax: float
            Maximum value of the Y axis.
            This is not Y parameter passed to xKeywords.
        yUnit: str
            Unit of the Y axis.
            Supported: '(K)', '(mK)'
            This is not Y parameter passed to xKeywords.
        colorMax: float
            Maximum value of the Y axis.
        colorUnit: str
            Unit of the Y axis.
            Supported: '(eV)', '(meV)', '($\\mu$eV) <- Without one backslash'
        """
        X = [element[0] for element in self.params]
        X = sorted(list(set(X)))

        Y = [element[1] for element in self.params]
        Y = sorted(list(set(Y)))

        Z = [element[2] for element in self.params]
        Z = sorted(list(set(Z)))

        gammaZPlot = []
        firstXPlot = []
        secondXPlot = []
        yPlot = []

        secondXCallback = None
        neighborGammasList = []
        neighborKeys = []
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

        # Set energy units and multipliers
        if yUnit == r"(K)":
            yMultiplier = 1
        elif yUnit == r"(mK)":
            yMultiplier = 1e3
        else:
            raise ValueError(f"Unknown yUnit: {yUnit}")

        # Set energy units and multipliers
        if colorUnit == r"(meV)":
            colorMultiplier = 1
        elif colorUnit == r"($\mu$eV)":
            colorMultiplier = 1e3
        elif colorUnit == r"(eV)":
            colorMultiplier = 1e-3
        else:
            raise ValueError(f"Unknown yUnit: {yUnit}")

        # Main plotting loop
        for nNeighborhood, gammaDict in enumerate(neighborGammasList):
            for key in self.symmetryKeys[gammaNeighorhoodLabels[nNeighborhood]]:
                for z in Z:
                    fig = plt.figure(figsize=(7, 5), dpi=400)
                    # Set up GridSpec (1 row, 1 column, with some spacing)
                    gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.75, bottom=0.25)
                    ax1 = fig.add_subplot(gs[0,0])

                    xPlot = []
                    xPlotFixedY = []
                    secondXPlot = []
                    yPlot = []
                    gammaColorPlot = []
                    for x in X:
                        for i in range(len(self.params)):
                            if int(self.params[i][2]) == z:
                                xPlot.append(self.params[i][0])
                                yPlot.append(self.params[i][1] * yMultiplier)
                                gammaColorPlot.append(np.abs(gammaDict[key][i]) * colorMultiplier)
                            if self.params[i][1] == self.params[0][1]:
                                "To get carrier densities only at temperature == 0K"
                                xPlotFixedY.append(self.params[i][0] - firstXShift)
                                secondXPlot.append(secondXCallback(self.fillingTotal[i]))

                    #Creating a grid for colormap
                    xUnique = np.unique(xPlot)
                    yUnique = np.unique(yPlot)
                    xGrid, yGrid = np.meshgrid(xUnique, yUnique)

                    points = np.array([xPlot, yPlot]).T
                    Gap_grid = griddata(
                        points, gammaColorPlot, (xGrid, yGrid), method="linear", fill_value=0
                    )

                    #Plotting a colormap
                    colormesh = ax1.pcolormesh(
                        xGrid,
                        yGrid,
                        Gap_grid,
                        cmap="inferno",
                        norm=PowerNorm(gamma=.8, vmin=0.0, vmax=colorMax if colorMax != np.inf else max(gammaColorPlot)),
                    )

                    #Setting labels
                    band, spin1, spin2, sublat, symmetry = key
                    ax1.set_xlabel(firstXLabel)
                    ax1.set_ylabel(rf"T {yUnit}")
                    ax1.set_xlim(right=firstXMax if firstXMax != np.inf else max(xPlot))
                    ax1.xaxis.set_major_locator(ticker.LinearLocator(5))
                    ax1.yaxis.set_major_locator(ticker.LinearLocator(4))

                    colorbar = fig.colorbar(colormesh, ax=ax1)
                    colorbar.set_label(
                        rf"{gammaLabelsCallbacks[nNeighborhood](sublat, symmetry, spin1, spin2)}" + colorUnit,
                    )

                    # Getting second X axis
                    plt.draw()
                    ax1_ticks = ax1.get_xticks()
                    # Plot secondary axis for occupation
                    tick_labels = np.interp(
                        ax1_ticks, xPlotFixedY, secondXPlot
                    )  # Interpolate the mapping
                    ax2 = ax1.secondary_xaxis("top")
                    ax2.set_xticks(ax1_ticks)  # Use the same positions as `ef_plot`
                    ax2.set_xticklabels(
                        [f"{val:.1f}" for val in tick_labels]
                    )  # Map `n_total_plot` as tick labels
                    ax2.set_xlabel(secondXLabel, labelpad = 16)

                    ax1.tick_params(axis="x", direction="out")
                    ax1.tick_params(axis="y", direction="out")
                    ax2.tick_params(axis="x", direction="out")

                    # Saving figure
                    plt.savefig(
                        f"../Plots/GammaCmap_{gammaNeighorhoodLabels[nNeighborhood]}_{z}_band{band}_spin1{spin1}_spin2{spin2}_lat{sublat}_{symmetry}.png"
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

        ax1.legend(title=r"$U$ (meV)", loc="best")
        ax1.set_xlabel(r"$\mu$ (meV)")
        ax1.set_ylabel(r"$n$ (10 \textsuperscript{-2})", labelpad=10)
        ax1.set_xlim(right=150)
        ax1.set_ylim(top=16)
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(30))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax1.grid(True, linestyle=':')
        plt.draw()

        ax1_ticks = ax1.get_yticks()
        tick_labels = []
        for tick in ax1_ticks:
            tick_labels.append(self.__calculateFillingPerCm2(tick / 100 * 12))
        # Plot secondary axis for occupation
        # ax2 = ax1.secondary_yaxis("right")
        # ax2.set_yticks(ax1_ticks)  # Use the same positions as `ef_plot`
        # ax2.set_yticklabels([f"{val:.1f}" for val in tick_labels])  # Map `n_total_plot` as tick labels
        # ax2.set_ylabel(r"$n$ (10 \textsuperscript{-14} cm\textsuperscript{-2})", labelpad=45, rotation=270)  # Customize units as needed


        plt.savefig(f"../Plots/FillingFermiTotal.png")
        plt.close()

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __initializeSymmetryKeys(self):
        # Nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin1 in range(1, 3):
                for spin2 in range(1, 3):
                    for sublat in range(1, self.layerCouplings + 1):
                        for symmetry in self.projector.getSymmetryNames():
                            self.symmetryKeys["nearest"].append(
                                (band, spin1, spin2, sublat, symmetry)
                            )

        # Next-to-nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin1 in range(1, 3):
                for spin2 in range(1, 3):
                    for sublat in range(1, self.sublattices + 1):
                        for symmetry in self.projector.getSymmetryNames():
                            self.symmetryKeys["next"].append(
                                (band, spin1, spin2, sublat, symmetry)
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
        self, sublat: int, symmetry: str, spin1: int, spin2: int
    ) -> str:
        return rf"$\Gamma_{{\alpha \overline{{\alpha}}}}^{{{symmetry}}}$"

    def __getNextNearestNeighborGammaLabel(
        self, sublat: int, symmetry: str, spin1: int, spin2: int
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
