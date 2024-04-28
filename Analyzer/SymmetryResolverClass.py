import numpy as np
from DataReaderClass import *

class SymmetryResolver(DataReader):
    
    def __init__(self, nNeighbors: int, runsPath: str, matchPattern: str):
        """
        Initializes SymmetryResolver object to callculate supeconducting gap symmetries.
        Arguments:
            self.nNeighbors - number of neighbors that are paired to each atomic site. 
                              For hexagonal lattice we have 3 nearest neighbors and 6 next nearest neighbors.
            self.runsPath, self.matchPatterr - used to initialize DataReader object. See DataReade documentation.
        """
        DataReader.__init__(self, runsPath, matchPattern)
        self.nNeighbors = nNeighbors
        self.symmetryGammaDict: dict = {}

    #RETHINK ___WavePairing construction to avoid code repetition
    def _SWavePairing(self, listOfGammas: list) -> np.float64:
        """
        Sets parameters for s-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 0
        M = 0
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _PWavePairing(self, listOfGammas: list):
        """
        Sets parameters for p-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 1
        M = 1
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _DWavePairing(self, listOfGammas: list):
        """
        Sets parameters for s-wave gap symmetry calculation and invokes function _CalculateSingleGamma()
        which implements general equation for gap symmetry.
        """
        p = 0
        M = 2
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _CalculateSingleGamma(self, listOfGammas: list, p: int, M: int) -> np.float64:
        """
        Sums all gammas with proper phase factors, implementing general equation for given symmetry, determined by M and p.
        !!! Assumes that all neigbors in list have the same phase offset to the previous one !!!
        """
        #TODO: Rethink whether for next nearest neighbors pairing the condition of equal phase offset is valid.
        symmetryGamma = 0
        for i in range(self.nNeighbors):
            currentPhase = 2*np.pi / self.nNeighbors * i
            phaseFactor = np.exp(-1j*M*currentPhase)
            symmetryGamma += listOfGammas[i]*phaseFactor

        return symmetryGamma*(1j)**p/self.nNeighbors

    def CalculateSymmetryGamma(self):
        """
        Fills dict of superconducting gap symmetries, based on data from DataReader.gamma
        """

        #Loop over all dict keys except for neighbours
        for spin in [1,2]:
            for sublat in [1,2]:
                for orbital in [1,2,3]:

                    #Loop over all gammas for given spin, sublat and orbital
                    for i in range(len(self.gamma[(1,1,1,1)])):
                        gammaToSymmetrize = [] #THis code could be simplified, maybe list comprehension?
                        for neighbor in range(1, self.nNeighbors + 1):
                            gammaToSymmetrize.append(self.gamma[(spin, neighbor, sublat, orbital)][i])
                        if i == 0:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')] = [self._SWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p')] = [self._PWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd')] = [self._DWavePairing(gammaToSymmetrize)]
                        else:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')].append(self._SWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p')].append(self._PWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd')].append(self._DWavePairing(gammaToSymmetrize))

