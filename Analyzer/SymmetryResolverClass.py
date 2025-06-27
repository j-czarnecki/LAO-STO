import numpy as np
from DataReaderClass import *
from collections import defaultdict
from Projectors.Projectors import *

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
        self.projector = ProjectorC3v() #This should be changed according to group based on which symetry is to be analyzed

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
            for spin in range(1, 3):
                for sublat in range(1, self.layerCouplings + 1):
                    symmetryKey = (band, spin, sublat)

                    for i in range(listLen):
                        gammaToSymmetrize = []
                        for orbital in range(1, 4):
                            # This code could be simplified, maybe list comprehension?
                            # Nearest neighbours pairing
                            for sublat in range(1, self.layerCouplings + 1):
                                for neighbor in range(1, self.nNeighbors + 1):
                                    gammaKey = (
                                        (spin, neighbor, sublat, orbital)
                                        if self.subbands == 0
                                        else (band, spin, neighbor, sublat, orbital)
                                    )
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
            for spin in range(1, 3):
                for sublat in range(1, self.sublattices + 1):
                    symmetryKey = (band, spin, sublat)
                    for i in range(listLen):
                        nnnGammaToSymmetrize = []
                        for orbital in range(1, 4):
                            # Next nearest neighbours pairing
                            if self.nNextNeighbors != 0:
                                for neighbor in range(
                                    self.nNeighbors + 1,
                                    self.nNeighbors + self.nNextNeighbors + 1,
                                ):
                                    gammaKey = (
                                        (spin, neighbor, sublat, orbital)
                                        if self.subbands == 0
                                        else (band, spin, neighbor, sublat, orbital)
                                    )
                                    nnnGammaToSymmetrize.append(self.gamma[gammaKey][i])
                        for i, sym in enumerate(symmetries):
                            self.nnnSymmetryGammaDict[(*symmetryKey, sym)].append(
                                symmetryCallbacks[i](nnnGammaToSymmetrize, "next")
                            )

    def calculateSingletTripletGammas(self) -> None:
        if not self.symmetryGammaDict:
            raise ValueError(
                "Dictionary self.symmetryGammaSingletTripletDict is empty. Please, calculate symmetries first."
            )
        if self.nNextNeighbors != 0 and not self.nnnSymmetryGammaDict:
            raise ValueError(
                "Dictionary self.nnnSymmetryGammaSingletTripletDict is empty. Please, calculate symmetries first."
            )

        for key in self.symmetryGammaDict:
            band, spin, sublat, orbital, symmetry = key
            spinOpposite = spin + (-1) ** (spin + 1)  # generates opposite spin
            spinStateSign = (-1) ** spin  # -1 for singlet state, +1 for triplet state
            keyOpposite = (band, spinOpposite, sublat, orbital, symmetry)

            self.symmetryGammaSingletTripletDict[key] = []
            for i in range(len(self.symmetryGammaDict[key])):
                value = 0.5 * (  # Check whether it should be 0.5 or 1/sqrt(2)
                    self.symmetryGammaDict[key][i]
                    + np.conj(self.symmetryGammaDict[keyOpposite][i])
                    * spinStateSign
                )
                self.symmetryGammaSingletTripletDict[key].append(value)

        if self.nNextNeighbors == 0:
            return

        for key in self.nnnSymmetryGammaDict:
            band, spin, sublat, orbital, symmetry = key
            spinOpposite = spin + (-1) ** (spin + 1)  # generates opposite spin
            spinStateSign = (-1) ** spin  # -1 for singlet state, +1 for triplet state
            keyOpposite = (band, spinOpposite, sublat, orbital, symmetry)

            self.nnnSymmetryGammaSingletTripletDict[key] = []
            for i in range(len(self.nnnSymmetryGammaDict[key])):
                value = 0.5 * (  # Check whether it should be 0.5 or 1/sqrt(2)
                    self.nnnSymmetryGammaDict[key][i]
                    + np.conj(self.nnnSymmetryGammaDict[keyOpposite][i]) * spinStateSign
                )
                self.nnnSymmetryGammaSingletTripletDict[key].append(value)

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
