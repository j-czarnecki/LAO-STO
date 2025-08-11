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