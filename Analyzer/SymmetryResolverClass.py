import numpy as np
from DataReaderClass import *

class SymmetryResolver(DataReader):
    
    def __init__(self, nNeighbors: int, nNextNeighbors: int, runsPath: str, matchPattern: str):
        """
        Initializes SymmetryResolver object to callculate supeconducting gap symmetries.
        Arguments:
            self.nNeighbors - number of neighbors that are paired to each atomic site. 
                              For hexagonal lattice we have 3 nearest neighbors.
            self.nNeighbors - number of next nearest neighbors that are paired to each atomic site. 
                              For hexagonal lattice we have 6 next nearest neighbors.
            self.runsPath, self.matchPatterr - used to initialize DataReader object. See DataReade documentation.
        """
        DataReader.__init__(self, runsPath, matchPattern)
        self.nNeighbors = nNeighbors
        self.nNextNeighbors = nNextNeighbors
        self.symmetryGammaDict: dict = {}
        self.nnnSymmetryGammaDict: dict = {}

    #RETHINK ___WavePairing construction to avoid code repetition
    def _ExtSWavePairing(self, listOfGammas: list) -> np.float64:
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
    
    def _CalculateSingleGamma(self, listOfGammas: list, p: int, M: int) -> np.float64:
        """
        Sums all gammas with proper phase factors, implementing general equation for given symmetry, determined by M and p.
        !!! Assumes that all neigbors in list have the same phase offset to the previous one !!!
        """
        symmetryGamma = 0
        nBonds = len(listOfGammas)
        for i in range(nBonds):
            currentPhase = 2*np.pi / nBonds * i
            phaseFactor = np.exp(-1j*M*currentPhase)
            symmetryGamma += listOfGammas[i]*phaseFactor

        return symmetryGamma*(1j)**p/nBonds

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
                        #Nearest neighbours pairing
                        for neighbor in range(1, self.nNeighbors + 1):
                            gammaToSymmetrize.append(self.gamma[(spin, neighbor, sublat, orbital)][i])
                        if i == 0:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')] = [self._ExtSWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p+')] = [self._PPlusWavePairing(gammaToSymmetrize)]
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd+')] = [self._DPlusWavePairing(gammaToSymmetrize)]
                        else:
                            self.symmetryGammaDict[(spin,sublat,orbital, 's')].append(self._ExtSWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'p+')].append(self._PPlusWavePairing(gammaToSymmetrize))
                            self.symmetryGammaDict[(spin,sublat,orbital, 'd+')].append(self._DPlusWavePairing(gammaToSymmetrize))


                        nnnGammaToSymmetrize = [] #THis code could be simplified, maybe list comprehension?
                        #Next nearest neighbours pairing
                        for neighbor in range(self.nNeighbors + 1, self.nNeighbors + self.nNextNeighbors + 1):
                            nnnGammaToSymmetrize.append(self.gamma[(spin, neighbor, sublat, orbital)][i])
                        if i == 0:
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 's')] = [self._ExtSWavePairing(nnnGammaToSymmetrize)]
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'p+')] = [self._PPlusWavePairing(nnnGammaToSymmetrize)]
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'd+')] = [self._DPlusWavePairing(nnnGammaToSymmetrize)]
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'f')] = [self._FWavePairing(nnnGammaToSymmetrize)]
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'd-')] = [self._DMinusWavePairing(nnnGammaToSymmetrize)]
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'p-')] = [self._PMinusWavePairing(nnnGammaToSymmetrize)]
                        else:
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 's')].append(self._ExtSWavePairing(nnnGammaToSymmetrize))
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'p+')].append(self._PPlusWavePairing(nnnGammaToSymmetrize))
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'd+')].append(self._DPlusWavePairing(nnnGammaToSymmetrize))
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'f')].append(self._FWavePairing(nnnGammaToSymmetrize))
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'd-')].append(self._DMinusWavePairing(nnnGammaToSymmetrize))
                            self.nnnSymmetryGammaDict[(spin,sublat,orbital, 'p-')].append(self._PMinusWavePairing(nnnGammaToSymmetrize))

                        

