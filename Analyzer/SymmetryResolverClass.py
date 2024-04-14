import numpy as np
from DataReaderClass import *

class SymmetryResolver(DataReader):
    
    def __init__(self, nNeighbors, runsPath, matchPattern):
        DataReader.__init__(self, runsPath, matchPattern)
        self.nNeighbors = nNeighbors
        self.symmetryGammaDict = {}

    #RETHINK ___WavePairing construction to avoid code repetition
    def _SWavePairing(self, listOfGammas: list):
        p = 0
        M = 0
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _PWavePairing(self, listOfGammas: list):
        p = 1
        M = 1
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _DWavePairing(self, listOfGammas: list):
        p = 0
        M = 2
        return self._CalculateSingleGamma(listOfGammas, p, M)

    def _CalculateSingleGamma(self, listOfGammas: list, p: int, M: int):
        symmetryGamma = 0
        for i in range(self.nNeighbors):
            currentPhase = 2*np.pi / self.nNeighbors * i
            phaseFactor = np.exp(-1j*M*currentPhase)
            symmetryGamma += listOfGammas[i]*phaseFactor

        return symmetryGamma*(1j)**p/self.nNeighbors

    def CalculateSymmetryGamma(self):

        #Loop over all dict keys except from neighbours
        for spin in [1,2]:
            for sublat in [1,2]:
                for orbital in [1,2,3]:

                    #Loop over all gammas for given spin, sublat and orbital
                    for i in range(len(self.gamma[(1,1,1,1)])):
                        gammaToSymmetrize = []
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

