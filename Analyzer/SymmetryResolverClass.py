import numpy as np
from DataReaderClass import *
from collections import defaultdict

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
            self.runsPath, self.matchPatterr - used to initialize DataReader object. See DataReade documentation.
        """
        DataReader.__init__(self, runsPath, matchPattern, sublattices, subbands)
        self.nNeighbors = nNeighbors
        self.nNextNeighbors = nNextNeighbors
        self.symmetryGammaDict: dict = defaultdict(list)
        self.symmetryGammaSingletTripletDict: dict = defaultdict(list)
        self.nnnSymmetryGammaDict: dict = defaultdict(list)
        self.nnnSymmetryGammaSingletTripletDict: dict = defaultdict(list)

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def CalculateSymmetryGamma(self):
        """
        Fills dict of superconducting gap symmetries, based on data from DataReader.gamma
        """

        symmetries = [r"A_1^{(1)}", r"A_1^{(2)}",
                                     r"A_2^{(1)}",
                                     r"B_1^{(1)}", r"B_1^{(2)}",
                                     r"B_2^{(1)}",
                                     r"E_1^{(1)}", r"E_1^{(2)}", r"E_1^{(3)}", r"E_1^{(4)}", r"E_1^{(5)}", r"E_1^{(6)}",
                                     r"E_2^{(1)}", r"E_2^{(2)}", r"E_2^{(3)}", r"E_2^{(4)}", r"E_2^{(5)}", r"E_2^{(6)}"]
        symmetryCallbacks = [
            #s-wave
            self.__A1Projection1,
            self.__A1Projection2,
            #XXX-wave
            self.__A2Projection,
            #f-wave
            self.__B1Projection1,
            self.__B1Projection2,
            #YYY-wave
            self.__B2Projection,
            #p-wave
            self.__E1Projection1,
            self.__E1Projection2,
            self.__E1Projection3,
            self.__E1Projection4,
            self.__E1Projection5,
            self.__E1Projection6,
            #d-wave
            self.__E2Projection1,
            self.__E2Projection2,
            self.__E2Projection3,
            self.__E2Projection4,
            self.__E2Projection5,
            self.__E2Projection6,
        ]

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
                            for neighbor in range(1, self.nNeighbors + 1):
                                gammaKey = (
                                    (spin, neighbor, sublat, orbital)
                                    if self.subbands == 0
                                    else (band, spin, neighbor, sublat, orbital)
                                )
                                gammaToSymmetrize.append(self.gamma[gammaKey][i])
                                # Generating next atom from second sublattice
                                # such that it obeys the order of neighbors defined on a regular hexagon
                                gammaKeySecondSublat = (
                                    (spin, self.__getFollowingNeighbor(neighbor), self.__getOppositeSublat(sublat), orbital)
                                    if self.subbands == 0
                                    else (band, spin, self.__getFollowingNeighbor(neighbor), self.__getOppositeSublat(sublat), orbital)
                                )
                                gammaToSymmetrize.append(self.gamma[gammaKeySecondSublat][i])
                        for i, sym in enumerate(symmetries):
                            self.symmetryGammaDict[(*symmetryKey, sym)].append(
                                symmetryCallbacks[i](gammaToSymmetrize)
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
                                symmetryCallbacks[i](nnnGammaToSymmetrize)
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

    def __getProj(self, listOfGammas: list, indecesPlus: list, indecesMinus: list) -> np.complex128:
        res = 0
        for i in indecesPlus:
            res += listOfGammas[i]
        for i in indecesMinus:
            res -= listOfGammas[i]
        return res/len(listOfGammas)

    # s-wave
    def __A1Projection1(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [2, 5, 7, 10, 12, 15]
        return self.__getProj(listOfGammas, indecesPlus, [])

    def __A1Projection2(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [0, 1, 3, 4, 6, 8, 9, 11, 13, 14, 16, 17]
        return self.__getProj(listOfGammas, indecesPlus, [])

    def __A2Projection(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [1, 4, 6, 9, 14, 17]
        indecesMinus = [0, 3, 8, 11, 13, 16]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    # f-wave
    def __B1Projection1(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [5, 7, 15]
        indecesMinus = [2, 10, 12]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __B1Projection2(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [1, 3, 9, 11, 13, 17]
        indecesMinus = [0, 4, 6, 8, 14, 16]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __B2Projection(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [0, 1, 8, 9, 16, 17]
        indecesMinus = [3, 4, 6, 11, 13, 14]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    # p-wave
    def __E1Projection1(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [4, 9]
        indecesMinus = [1, 6]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __E1Projection2(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [5, 10]
        indecesMinus = [2, 7]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __E1Projection3(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [0, 11]
        indecesMinus = [3, 8]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __E1Projection4(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [2, 15]
        indecesMinus = [5, 12]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E1Projection5(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [3, 16]
        indecesMinus = [0, 13]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E1Projection6(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [4, 17]
        indecesMinus = [1, 14]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    # d-wave
    def __E2Projection1(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [6, 9]
        indecesMinus = [1, 4]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)

    def __E2Projection2(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [7, 10]
        indecesMinus = [2, 5]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E2Projection3(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [8, 11]
        indecesMinus = [0, 3]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E2Projection4(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [12, 15]
        indecesMinus = [2, 5]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E2Projection5(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [13, 16]
        indecesMinus = [0, 3]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)


    def __E2Projection6(self, listOfGammas: list) -> np.complex128:
        indecesPlus = [14, 17]
        indecesMinus = [1, 4]
        return self.__getProj(listOfGammas, indecesPlus, indecesMinus)



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
