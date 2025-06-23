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
    for irrep in self.irrepsTuple:
      print(f"IR: {irrep}")
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
            display(f"indecesPlus = {oneIndeces}")
            display(f"indecesMinus = {minusOneIndeces}")
        else:
          print("No matching eigenvalue")
      print("------------")
      print("\n")

  def __instantProjectionCheck(self, projection: sp.Matrix, irrep: str):
    """ Checks corectness of single projection operator """
    if ((projection**2 - projection) == np.zeros(projection.shape[0])):
      display(projection**2 - projection)
      raise Exception(f"Projection for irrep {irrep} operator is not indempotent")