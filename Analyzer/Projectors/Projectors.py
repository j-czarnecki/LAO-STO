import numpy as np
from abc import ABC, abstractmethod
from typing import Callable

class Projector(ABC):
  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  @abstractmethod
  def getCallbacks(self) -> list[Callable[[list[np.complex128]], np.complex128]]:
      pass

  @abstractmethod
  def getSymmetryNames(self) -> tuple[str,...]:
      pass

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Protected methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  def _getProj(self, listOfGammas: list[np.complex128], indecesPlus: list[int], indecesMinus: list[int]) -> np.complex128:
      res: np.complex128 = np.complex128(0)
      for i in indecesPlus:
          res += listOfGammas[i]
      for i in indecesMinus:
          res -= listOfGammas[i]
      return res/len(listOfGammas)

class ProjectorC6v(Projector):

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """
  def getSymmetryNames(self) -> tuple[str, ...]:
    symmetries = (r"A_1^{(1)}", r"A_1^{(2)}",
                  r"A_2^{(1)}",
                  r"B_1^{(1)}", r"B_1^{(2)}",
                  r"B_2^{(1)}",
                  r"E_1^{(1)}", r"E_1^{(2)}", r"E_1^{(3)}", r"E_1^{(4)}", r"E_1^{(5)}", r"E_1^{(6)}",
                  r"E_2^{(1)}", r"E_2^{(2)}", r"E_2^{(3)}", r"E_2^{(4)}", r"E_2^{(5)}", r"E_2^{(6)}")
    return symmetries

  
  def getCallbacks(self) -> list[Callable[[list[np.complex128]], np.complex128]]:
      symmetryCallbacks = [
          #s-wave
          self.A1Projection1,
          self.A1Projection2,
          #XXX-wave
          self.A2Projection,
          #f-wave
          self.B1Projection1,
          self.B1Projection2,
          #YYY-wave
          self.B2Projection,
          #p-wave
          self.E1Projection1,
          self.E1Projection2,
          self.E1Projection3,
          self.E1Projection4,
          self.E1Projection5,
          self.E1Projection6,
          #d-wave
          self.E2Projection1,
          self.E2Projection2,
          self.E2Projection3,
          self.E2Projection4,
          self.E2Projection5,
          self.E2Projection6,
      ]
      return symmetryCallbacks

  # s-wave
  def A1Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [2, 5, 7, 10, 12, 15]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A1Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [0, 1, 3, 4, 6, 8, 9, 11, 13, 14, 16, 17]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A2Projection(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [1, 4, 6, 9, 14, 17]
      indecesMinus = [0, 3, 8, 11, 13, 16]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  # f-wave
  def B1Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [5, 7, 15]
      indecesMinus = [2, 10, 12]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def B1Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [1, 3, 9, 11, 13, 17]
      indecesMinus = [0, 4, 6, 8, 14, 16]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def B2Projection(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [0, 1, 8, 9, 16, 17]
      indecesMinus = [3, 4, 6, 11, 13, 14]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  # p-wave
  def E1Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [4, 9]
      indecesMinus = [1, 6]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def E1Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [5, 10]
      indecesMinus = [2, 7]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def E1Projection3(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [0, 11]
      indecesMinus = [3, 8]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def E1Projection4(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [2, 15]
      indecesMinus = [5, 12]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E1Projection5(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [3, 16]
      indecesMinus = [0, 13]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E1Projection6(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [4, 17]
      indecesMinus = [1, 14]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  # d-wave
  def E2Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [6, 9]
      indecesMinus = [1, 4]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def E2Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [7, 10]
      indecesMinus = [2, 5]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E2Projection3(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [8, 11]
      indecesMinus = [0, 3]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E2Projection4(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [12, 15]
      indecesMinus = [2, 5]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E2Projection5(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [13, 16]
      indecesMinus = [0, 3]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def E2Projection6(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [14, 17]
      indecesMinus = [1, 4]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


class ProjectorC3v(Projector):

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """
  def getSymmetryNames(self) -> tuple[str, ...]:
    symmetries = (r"A_1^{(1)}", r"A_1^{(2)}", r"A_1^{(3)}", r"A_1^{(4)}",
                  r"A_2^{(1)}", r"A_2^{(2)}"
                  r"E^{(1)}", r"E^{(2)}", r"E^{(3)}", r"E^{(4)}", r"E^{(5)}", r"E^{(6)}",
                  r"E^{(7)}", r"E^{(8)}", r"E^{(9)}", r"E^{(10)}", r"E^{(11)}", r"E^{(12)}")
    return symmetries

  
  def getCallbacks(self) -> list[Callable[[list[np.complex128]], np.complex128]]:
      symmetryCallbacks = [
        self.A1Projection1,
        self.A1Projection2,
        self.A1Projection3,
        self.A1Projection4,
        self.A2Projection1,
        self.A2Projection2,
        self.EProjection1,
        self.EProjection2,
        self.EProjection3,
        self.EProjection4,
        self.EProjection5,
        self.EProjection6,
        self.EProjection7,
        self.EProjection8,
        self.EProjection9,
        self.EProjection10,
        self.EProjection11,
        self.EProjection12
      ]
      return symmetryCallbacks

  # s-wave
  def A1Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [2, 10, 12]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A1Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [5, 7, 15]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A1Projection3(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [0, 4, 6, 8, 14, 16]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A1Projection4(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [1, 3, 9, 11, 13, 17]
      return self._getProj(listOfGammas, indecesPlus, [])

  def A2Projection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [0, 8, 16]
      indecesMinus = [4, 6, 14]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def A2Projection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [1, 9, 17]
      indecesMinus = [3, 11, 13]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection1(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [6]
      indecesMinus = [4]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection2(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [7]
      indecesMinus = [5]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)


  def EProjection3(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [8]
      indecesMinus = [0]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection4(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [9]
      indecesMinus = [1]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection5(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [10]
      indecesMinus = [2]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection6(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [11]
      indecesMinus = [3]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection7(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [12]
      indecesMinus = [2]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection8(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [13]
      indecesMinus = [3]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection9(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [14]
      indecesMinus = [4]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection10(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [15]
      indecesMinus = [5]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection11(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [16]
      indecesMinus = [0]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)

  def EProjection12(self, listOfGammas: list[np.complex128]) -> np.complex128:
      indecesPlus = [17]
      indecesMinus = [1]
      return self._getProj(listOfGammas, indecesPlus, indecesMinus)




