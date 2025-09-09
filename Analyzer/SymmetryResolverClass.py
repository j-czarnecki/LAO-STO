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

import numpy as np
from DataReaderClass import *
from collections import defaultdict
from Projectors.Projectors import *
import logging

logger = logging.getLogger(__name__)


class SymmetryResolver(DataReader):

    def __init__(
        self,
        nNeighbors: int,
        nNextNeighbors: int,
        runsPath: str,
        matchPattern: str,
        sublattices: int,
        subbands: int,
    ):
        """
        Initializes SymmetryResolver object to callculate supeconducting gap symmetries.
        Arguments:
            self.nNeighbors - number of neighbors that are paired to each atomic site.
                              For hexagonal lattice we have 3 nearest neighbors.
            self.nNeighbors - number of next nearest neighbors that are paired to each atomic site.
                              For hexagonal lattice we have 6 next nearest neighbors.
            self.runsPath, self.matchPattern - used to initialize DataReader object. See DataReade documentation.
        """
        DataReader.__init__(self, runsPath, matchPattern, sublattices, subbands)
        self.nNeighbors = nNeighbors
        self.nNextNeighbors = nNextNeighbors
        self.symmetryGammaDict: dict = defaultdict(list)
        self.symmetryGammaSingletTripletDict: dict = defaultdict(list)
        self.nnnSymmetryGammaDict: dict = defaultdict(list)
        self.nnnSymmetryGammaSingletTripletDict: dict = defaultdict(list)
        self.projector = ProjectorC6v() #This should be changed according to group based on which symetry is to be analyzed

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def CalculateSymmetryGamma(self):
        """
        Fills dict of superconducting gap symmetries, based on data from DataReader.gamma
        """

        symmetries = self.projector.getSymmetryNames()
        symmetryCallbacks = self.projector.getCallbacks()

        # Returns length of lists stored in dictionary
        listLen = len(next(iter(self.gamma.values())))

        # Loop over all dict keys except for neighbours
        # Nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin1 in range(1, 3):
                for spin2 in range(1, 3):
                    for sublat in range(1, self.layerCouplings + 1):
                        symmetryKey = (band, spin1, spin2, sublat)

                        for i in range(listLen):
                            gammaToSymmetrize = []
                            for orbital in range(1, 4):
                                # This code could be simplified, maybe list comprehension?
                                # Nearest neighbours pairing
                                for sublat in range(1, self.layerCouplings + 1):
                                    for neighbor in range(1, self.nNeighbors + 1):
                                        gammaKey = (band, spin1, spin2, neighbor, sublat, orbital)
                                        gammaToSymmetrize.append(self.gamma[gammaKey][i])
                                        # Generating next atom from second sublattice
                                        # such that it obeys the order of neighbors defined on a regular hexagon
                                        # gammaKeySecondSublat = (
                                        #     (spin, self.__getFollowingNeighbor(neighbor), self.__getOppositeSublat(sublat), orbital)
                                        #     if self.subbands == 0
                                        #     else (band, spin, self.__getFollowingNeighbor(neighbor), self.__getOppositeSublat(sublat), orbital)
                                        # )
                                        #gammaToSymmetrize.append(self.gamma[gammaKeySecondSublat][i])
                            for i, sym in enumerate(symmetries):
                                self.symmetryGammaDict[(*symmetryKey, sym)].append(
                                    symmetryCallbacks[i](gammaToSymmetrize, "nearest")
                                )
        if self.nNextNeighbors == 0:
            return

        # Next-to-nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin1 in range(1, 3):
                for spin2 in range(1, 3):
                    for sublat in range(1, self.sublattices + 1):
                        symmetryKey = (band, spin1, spin2, sublat)
                        for i in range(listLen):
                            nnnGammaToSymmetrize = []
                            for orbital in range(1, 4):
                                # Next nearest neighbours pairing
                                if self.nNextNeighbors != 0:
                                    for neighbor in range(
                                        self.nNeighbors + 1,
                                        self.nNeighbors + self.nNextNeighbors + 1,
                                    ):
                                        gammaKey = (band, spin1, spin2, neighbor, sublat, orbital)
                                        nnnGammaToSymmetrize.append(self.gamma[gammaKey][i])
                            for i, sym in enumerate(symmetries):
                                self.nnnSymmetryGammaDict[(*symmetryKey, sym)].append(
                                    symmetryCallbacks[i](nnnGammaToSymmetrize, "next")
                                )

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def __getOppositeSublat(self, sublat: int) -> int:
        """
        Returns opposite sublattice with which the coupling should be calculated
        """
        return sublat + (-1)**(sublat + 1)

    def __getFollowingNeighbor(self, neighbor: int) -> int:
        """
        Returns neighbors belonging to opposite lattice that should be takes as next
        in order determined by placing all the atoms on a regular hexagon for symetry analysis.
        """
        return (neighbor + 1) % self.nNeighbors + 1
    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """
