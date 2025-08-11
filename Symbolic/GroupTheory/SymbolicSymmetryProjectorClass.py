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

import sympy as sp
import numpy as np
from IPython.display import display


class SymbolicSymmetryProjectorClass:

  def __init__(self, irrepsTuple: tuple[str, ...], conjugacyClassesTuple: tuple[str, ...], characterTableDict: dict[str, dict[str, sp.Expr]], representationMatrixSize: int):
    """
    Initializes ScGapSymmetryProjector object to callculate supeconducting gap symmetries for different submodules.
    Arguments:
      irrepsTuple: tuple[str] - tuple of irreps in the analyzed point group.
      conjugacyClassesTuple: tuple[str] - tuple of conjugacy classes in the analyzed point group.
      characterTableDict: dict[str, dict[str, float]] - dictionary of character table, with characters
                                                        for each irrep (first key) and conjugacy class (second key).
      representationMatrixSize: int - size of representation matrix in given submodule.
    Initializes:
      ...
    """
    self.irrepsTuple: tuple[str, ...] = irrepsTuple
    self.conjugacyClassesTuple: tuple[str, ...] = conjugacyClassesTuple
    self.characterTableDict: dict[str, dict[str, sp.Expr]] = characterTableDict
    self.representationMatrixSize: int = representationMatrixSize
    self.g = sum(characterTableDict[irrep][conjugacyClassesTuple[0]]**2 for irrep in irrepsTuple)

  def getProjectionOperators(self, operationsDict: dict[str, list[sp.Matrix]]) -> dict[str, sp.Matrix]:
    """
    Sets symmetry operations dictionary.
    Arguments:
      operationsDict: dict[str, list[float]] - dictionary of symmetry operations in each conjugacy class.
                                             The key is the name of conjugacy class and list contains
                                             all symmetry operations in that conjugacy class.
    """
    projections = {}
    for irrep in self.irrepsTuple:
      dGamma = self.characterTableDict[irrep][self.conjugacyClassesTuple[0]]
      projections[irrep] = sp.zeros(self.representationMatrixSize, self.representationMatrixSize)
      for conjugacyClass in self.conjugacyClassesTuple:
        for operation in operationsDict[conjugacyClass]:
          projections[irrep] += dGamma / self.g * self.characterTableDict[irrep][conjugacyClass] * operation

      self.__instantProjectionCheck(projections[irrep], irrep)

    return projections

  def getMultiplicities(self, operationsDict: dict[str, list[sp.Matrix]]) -> dict[str, int]:
    multiplicitiesDict = {}
    for irrep in self.irrepsTuple:
      multiplicitiesDict[irrep] = 0
      for conjugacyClass in self.conjugacyClassesTuple:
        for operation in operationsDict[conjugacyClass]:
          multiplicitiesDict[irrep] += self.characterTableDict[irrep][conjugacyClass] / self.g * operation.trace()

    return multiplicitiesDict

  def getDiagonalizedProjections(self, projectionsDict: dict[str, sp.Matrix]) -> dict[str, list[tuple[int, int, list[sp.Matrix]]]]:
    """
    Diagonalizes projection operators and return tuples of (eigenvalue, multiplicity, eigenvectors).
    Arguments:
      projectionsDict: dict[str, sp.Matrix] - dictionary of projection operators for each irrep (key).
    """
    diagonalizedProjections = {}
    for irrep in self.irrepsTuple:
      diagonalizedProjections[irrep] = projectionsDict[irrep].eigenvects()

    return diagonalizedProjections

  def displayProjectionsMetadata(self, diagonalizedProjections: dict[str, list[tuple[int, int, list[sp.Matrix]]]]):
    totalOutputStr: str = ""
    for irrep in self.irrepsTuple:
      print(f"IR: {irrep}")
      degeneracyCount = 0
      #Only interested in eigenvalues +-1
      for v in diagonalizedProjections[irrep]:
        if np.abs(v[0]) == 1:
          # Print indeces of +/- ones for easy python implementation
          for i in range(len(v[2])):
            oneIndeces = []
            minusOneIndeces = []
            for j in range(len(v[2][i])):
              if v[2][i][j] == 1:
                oneIndeces.append(j)
              elif v[2][i][j] == -1:
                minusOneIndeces.append(j)
            #Weird notation to enable easy copying to projectionIndeces dict
            totalOutputStr += f'r"{irrep[:-1]}^{{({i + 1})}}$" : {{"plus" : {oneIndeces}, "minus" : {minusOneIndeces} }},\n'
        else:
          print("No matching eigenvalue")
      print("------------")
      print("\n")
    print(totalOutputStr)

  def __instantProjectionCheck(self, projection: sp.Matrix, irrep: str):
    """ Checks corectness of single projection operator """
    if ((projection**2 - projection) == np.zeros(projection.shape[0])):
      display(projection**2 - projection)
      raise Exception(f"Projection for irrep {irrep} operator is not indempotent")