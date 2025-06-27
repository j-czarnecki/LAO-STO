import numpy as np
import pandas as pd
import os
import sys
import fortranformat as ff

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '/home', 'czarnecki', 'LAO-STO', 'Analyzer')))

from Projectors.Projectors import *

class OutputMocker:
  def __init__(self, outputPath: str, nOrbs: int = 3, nBands: int = 1, nSublats: int = 2):
    self.outputPath: str = outputPath
    self.nOrbs: int = nOrbs
    self.nBands: int = nBands
    self.nSublats: int = nSublats
    self.nNeighbors: int = 3
    self.nNextNeighbors: int = 6
    self.nSpins: int = 2
    self.layerCouplings = 2 * (self.nSublats - 1)
    self.projector = ProjectorC6v()

  def mockChargeOutput(self) -> None:
    mockChargeValue = 1.
    fortFormat = ff.FortranRecordWriter("(4I5, 1E15.5)")
    with open(os.path.join(self.outputPath, "OutputData", "Charge_dens_final.dat"), "w") as f:
      print(" #band spin lattice orbital Charge", file=f)
      for band in range(self.nBands):
        for spin in range(self.nSpins):
          for sublat in range(self.nSublats):
            for orb in range(self.nOrbs):
              line = fortFormat.write([band + 1, spin + 1, sublat + 1, orb + 1, mockChargeValue])
              print(line, file=f)

  def mockGammaOutput(self, gammaAmplitude: np.complex128, symmetries: tuple[str, ...]) -> None:

    fortFormat = ff.FortranRecordWriter("(5I5, 2E15.5)")
    with open(os.path.join(self.outputPath, "OutputData", "Gamma_SC_final.dat"), "w") as f:
      print(" #band spin neighbour lattice orbital Re(Gamma) Im(Gamma)", file=f)
      for band in range(self.nBands):
        for spin in range(self.nSpins):
          spinSign = (-1)**(spin) #Assuming spin singlet
          symGamma = self.__createFlatSymmetryGamma(symmetries) * gammaAmplitude * spinSign
          for neigh in range(self.nNeighbors + self.nNextNeighbors):
            maxLat = self.layerCouplings if neigh <= self.nNeighbors else self.nSublats
            for lat in range(maxLat):
              for orb in range(self.nOrbs):
                if neigh < self.nNeighbors:
                  gammaIdx = self.__getGammaIdx(orb, lat, neigh)
                  line = fortFormat.write([band + 1, spin + 1, neigh + 1, lat + 1, orb + 1, symGamma[gammaIdx].real, symGamma[gammaIdx].imag])
                else:
                  line = fortFormat.write([band + 1, spin + 1, neigh + 1, lat + 1, orb + 1, 0, 0])
                print(line, file=f)
            print(" ", file=f)
            print(" ", file=f)

  def __createFlatSymmetryGamma(self, symmetries: tuple[str, ...]) -> np.ndarray:
    projectionIndecesNearestDict = self.projector.getProjectionIndeces()["nearest"]
    gammaFlat = np.zeros((self.nOrbs * self.nSublats * self.nNeighbors ), dtype=np.complex128)
    for symmetry in symmetries:
      for idx in projectionIndecesNearestDict[symmetry]["plus"]:
        gammaFlat[idx] += 1
      for idx in projectionIndecesNearestDict[symmetry]["minus"]:
        gammaFlat[idx] -= 1

    return gammaFlat

  def __getGammaIdx(self, orb, lat, neigh):
    return orb * self.nSublats * self.nNeighbors + lat * self.nNeighbors + neigh



def main():
  mocker = OutputMocker("/home/czarnecki/LAO-STO/")
  mocker.mockChargeOutput()
  mocker.mockGammaOutput(np.complex128(0.1), (r"$A_1^{(1)}$",))

if __name__ == "__main__":
  main()