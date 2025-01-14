import numpy as np
from DataReaderClass import *


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
        self.symmetryGammaDict: dict = {}
        self.nnnSymmetryGammaDict: dict = {}

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Interface methods ----------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def CalculateSymmetryGamma(self):
        """
        Fills dict of superconducting gap symmetries, based on data from DataReader.gamma
        """

        # Loop over all dict keys except for neighbours
        # Nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.layerCouplings + 1):
                    for orbital in range(1, 4):
                        symmetryKey = (band, spin, sublat, orbital)
                        # Returns length of lists stored in dictionary
                        listLen = len(next(iter(self.gamma.values())))
                        # Loop over all gammas for given spin, sublat and orbital
                        for i in range(listLen):
                            # This code could be simplified, maybe list comprehension?
                            gammaToSymmetrize = []
                            # Nearest neighbours pairing
                            for neighbor in range(1, self.nNeighbors + 1):
                                gammaKey = (
                                    (spin, neighbor, sublat, orbital)
                                    if self.subbands == 0
                                    else (band, spin, neighbor, sublat, orbital)
                                )
                                gammaToSymmetrize.append(self.gamma[gammaKey][i])
                            if i == 0:
                                self.symmetryGammaDict[(*symmetryKey, "s")] = [
                                    self._ExtSWavePairing(gammaToSymmetrize)
                                ]
                                self.symmetryGammaDict[(*symmetryKey, "p+")] = [
                                    self._PPlusWavePairing(gammaToSymmetrize)
                                ]
                                self.symmetryGammaDict[(*symmetryKey, "d+")] = [
                                    self._DPlusWavePairing(gammaToSymmetrize)
                                ]
                            else:
                                self.symmetryGammaDict[(*symmetryKey, "s")].append(
                                    self._ExtSWavePairing(gammaToSymmetrize)
                                )
                                self.symmetryGammaDict[(*symmetryKey, "p+")].append(
                                    self._PPlusWavePairing(gammaToSymmetrize)
                                )
                                self.symmetryGammaDict[(*symmetryKey, "d+")].append(
                                    self._DPlusWavePairing(gammaToSymmetrize)
                                )
        # Next-to-nearest neighbors
        for band in range(1, max(1, self.subbands) + 1):
            for spin in range(1, 3):
                for sublat in range(1, self.sublattices + 1):
                    for orbital in range(1, 4):
                        symmetryKey = (band, spin, sublat, orbital)
                        # THis code could be simplified, maybe list comprehension?
                        listLen = len(next(iter(self.gamma.values())))
                        for i in range(listLen):
                            nnnGammaToSymmetrize = []
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
                                if i == 0:
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "s")] = [
                                        self._ExtSWavePairing(nnnGammaToSymmetrize)
                                    ]
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "p+")] = [
                                        self._PPlusWavePairing(nnnGammaToSymmetrize)
                                    ]
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "d+")] = [
                                        self._DPlusWavePairing(nnnGammaToSymmetrize)
                                    ]
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "f")] = [
                                        self._FWavePairing(nnnGammaToSymmetrize)
                                    ]
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "d-")] = [
                                        self._DMinusWavePairing(nnnGammaToSymmetrize)
                                    ]
                                    self.nnnSymmetryGammaDict[(*symmetryKey, "p-")] = [
                                        self._PMinusWavePairing(nnnGammaToSymmetrize)
                                    ]
                                else:
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "s")
                                    ].append(
                                        self._ExtSWavePairing(nnnGammaToSymmetrize)
                                    )
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "p+")
                                    ].append(
                                        self._PPlusWavePairing(nnnGammaToSymmetrize)
                                    )
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "d+")
                                    ].append(
                                        self._DPlusWavePairing(nnnGammaToSymmetrize)
                                    )
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "f")
                                    ].append(self._FWavePairing(nnnGammaToSymmetrize))
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "d-")
                                    ].append(
                                        self._DMinusWavePairing(nnnGammaToSymmetrize)
                                    )
                                    self.nnnSymmetryGammaDict[
                                        (*symmetryKey, "p-")
                                    ].append(
                                        self._PMinusWavePairing(nnnGammaToSymmetrize)
                                    )

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Private methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """

    def _ExtSWavePairing(self, listOfGammas: list):
        """
        Sets parameters for s-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 0
        M = 0
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _PPlusWavePairing(self, listOfGammas: list):
        """
        Sets parameters for p+ip-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 1
        M = 1
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _DPlusWavePairing(self, listOfGammas: list):
        """
        Sets parameters for d+id-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 0
        M = 2
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _FWavePairing(self, listOfGammas: list):
        """
        Sets parameters for f-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 1
        M = 3
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _DMinusWavePairing(self, listOfGammas: list):
        """
        Sets parameters for d-id-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 0
        M = 4
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _PMinusWavePairing(self, listOfGammas: list):
        """
        Sets parameters for d-id-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 1
        M = 5
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _CalculateSingleGamma(self, listOfGammas: list, p: int, M: int):
        """
        Sums all gammas with proper phase factors, implementing general equation for given symmetry, determined by M and p.
        !!! Assumes that all neigbors in list have the same phase offset to the previous one !!!
        """
        symmetryGamma = 0
        nBonds = len(listOfGammas)
        for i in range(nBonds):
            currentPhase = 2 * np.pi / nBonds * i
            phaseFactor = np.exp(-1j * M * currentPhase)
            symmetryGamma += listOfGammas[i] * phaseFactor

        return symmetryGamma * (1j) ** p / nBonds

    """ ---------------------------------------------------------------------------------- """
    """ ---------------------------- Special methods ------------------------------------- """
    """ ---------------------------------------------------------------------------------- """
