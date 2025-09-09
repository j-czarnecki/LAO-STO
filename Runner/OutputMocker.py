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
    """
    Constructor for OutputMocker class.
    Parameters:
    outputPath (str): Path to output directory.
    nOrbs (int): Number of orbitals.
    nBands (int): Number of bands.
    nSublats (int): Number of sublattices.
    """
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
    """
    Produces an artificial Charge_dens_final.dat file in a directory specified by self.outputPath.
    """
    mockChargeValue = 1. #Value is not relevant unless Hubbard interaction is turned on.
    fortFormat = ff.FortranRecordWriter("(4I5, 1E15.5)")
    with open(os.path.join(self.outputPath, "OutputData", "Charge_dens_final.dat"), "w") as f:
      print(" #band spin lattice orbital Charge", file=f)
      for band in range(self.nBands):
        for spin in range(self.nSpins):
          for sublat in range(self.nSublats):
            for orb in range(self.nOrbs):
              line = fortFormat.write([band + 1, spin + 1, sublat + 1, orb + 1, mockChargeValue])
              print(line, file=f)

  def mockGammaOutput(self, gammaAmplitudeDict: dict[str, np.complex128], symmetriesWeightsDict: dict[str, dict[str, float]]) -> None:
    """
    Produces an artificial Gamma_SC_final.dat file in a directory specified by self.outputPath.
    Resultng Gamma is a linear combination of symmetries specified in symmetriesWeightsDict.
    Magnitude of the Gamma is specified by gammaAmplitudeDict.

    Parameters:
    gammaAmplitudeDict (dict[str, np.complex128]): Dictionary specifying amplitudes of symmetries for nearest and next-nearest neighbors.
    symmetriesWeightsDict (dict[str, dict[str, float]]): Dictionary specifying weights of symmetries for nearest and next-nearest neighbors.
    """

    fortFormat = ff.FortranRecordWriter("(6I5, 2E15.5)")
    with open(os.path.join(self.outputPath, "OutputData", "Gamma_SC_final.dat"), "w") as f:
      print(" #band spin neighbour lattice orbital Re(Gamma) Im(Gamma)", file=f)
      for band in range(self.nBands):
        for spin1 in range(self.nSpins):
          for spin2 in range(self.nSpins):
            if spin1 == spin2:
              continue
            spinSign = (-1)**(spin1) #Assuming spin singlet
            symGammaNearest = self.__createFlatSymmetryGamma(symmetriesWeightsDict["nearest"], "nearest") * gammaAmplitudeDict["nearest"] * spinSign
            symGammaNext = self.__createFlatSymmetryGamma(symmetriesWeightsDict["next"], "next") * gammaAmplitudeDict["next"] * spinSign
            for neigh in range(self.nNeighbors + self.nNextNeighbors):
              maxLat = self.layerCouplings if neigh <= self.nNeighbors else self.nSublats
              for lat in range(maxLat):
                for orb in range(self.nOrbs):
                  #Nearest neighbors
                  if neigh < self.nNeighbors:
                    gammaIdx = self.__getGammaIdx(orb, lat, neigh, "nearest")
                    line = fortFormat.write([band + 1, spin1 + 1, spin2 + 1, neigh + 1, lat + 1, orb + 1, symGammaNearest[gammaIdx].real, symGammaNearest[gammaIdx].imag])
                  #Next-nearest neighbors
                  else:
                    gammaIdx = self.__getGammaIdx(orb, lat, neigh - self.nNeighbors, "next") #Enumerate next-nearest neighbors from 0 to 6 for index getter
                    line = fortFormat.write([band + 1, spin1 + 1, spin2 + 1, neigh + 1, lat + 1, orb + 1, symGammaNext[gammaIdx].real, symGammaNext[gammaIdx].imag])
                  print(line, file=f)
              print(" ", file=f)
              print(" ", file=f)

  def __createFlatSymmetryGamma(self, symmetriesWeightsDict: dict[str, float], neighborsType: str = "nearest") -> np.ndarray:
    """
    Creates a flat array of gamma, taking into account weights of symmetries specified in symmetriesWeightsDict
    Constraint: Only two lattices currently supported
    """
    projectionIndecesDict = self.projector.getProjectionIndeces()[neighborsType]
    if neighborsType == "next":
      maxNeighbors = 6
      nCouplings = 1 #Because we do not take second sublat to irreps
    elif neighborsType == "nearest":
      maxNeighbors = 3
      nCouplings = self.layerCouplings
    else:
      raise ValueError("Unknown neighbors type")

    gammaFlat = np.zeros((self.nOrbs * nCouplings * maxNeighbors ), dtype=np.complex128)
    for symmetry in symmetriesWeightsDict.keys():
      for idx in projectionIndecesDict[symmetry]["plus"]:
        gammaFlat[idx] += 1 * symmetriesWeightsDict[symmetry]
      for idx in projectionIndecesDict[symmetry]["minus"]:
        gammaFlat[idx] -= 1 * symmetriesWeightsDict[symmetry]

    return gammaFlat

  def __getGammaIdx(self, orb: int, lat: int, neigh: int, neighborsType: str = "nearest"):
    """
    Calculates index of flattened gamma array for given orb, lat, neigh
    Constraint: Only two lattices currently supported
    """
    maxNeighbors = 1
    nCouplings = 1
    if neighborsType == "next":
      maxNeighbors = 6
      nCouplings = 1 #Because we do not take second sublat to irreps
      return orb * nCouplings * maxNeighbors + neigh
    elif neighborsType == "nearest":
      maxNeighbors = 3
      nCouplings = self.layerCouplings
      return orb * nCouplings * maxNeighbors + lat * maxNeighbors + neigh
    else:
      raise ValueError("Unknown neighbors type")



def main():
  mocker = OutputMocker("/home/czarnecki/LAO-STO/")
  mocker.mockChargeOutput()
  symmetriesWeightsDict = {"nearest": {
                              r"$A_1^{(1)}$": 0.1,
                              r"$A_1^{(2)}$": 0.4,
                            },
                          "next": {
                              r"$A_1^{(1)}$": 0.2,
                              r"$A_1^{(2)}$": 0.8,
                            },
                          }
  gammaAmplitudesDict = {"nearest": np.complex128(0.1),
                         "next": np.complex128(0.1)}
  mocker.mockGammaOutput(gammaAmplitudesDict, symmetriesWeightsDict)

if __name__ == "__main__":
  main()