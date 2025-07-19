import numpy as np
from abc import ABC, abstractmethod
from typing import Callable

class Projector(ABC):
  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  @abstractmethod
  def getCallbacks(self) -> list[Callable[[list[np.complex128], str], np.complex128]]:
      pass

  @abstractmethod
  def getSymmetryNames(self) -> tuple[str,...]:
      pass

  @abstractmethod
  def getProjectionIndeces(self) -> dict[str, dict[str, dict[str,list[int]]]]:
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

  def __init__(self):
    self.projectionIndeces: dict[str, dict[str, dict[str,list[int]]]] = {
        "nearest": {
            r"$A_1^{(1)}$" : {"plus" : [1, 4, 8, 11, 12, 15], "minus" : [] },
            r"$A_1^{(2)}$" : {"plus" : [0, 2, 3, 5, 6, 7, 9, 10, 13, 14, 16, 17], "minus" : [] },
            r"$A_2^{(1)}$" : {"plus" : [0, 3, 7, 10, 14, 17], "minus" : [2, 5, 6, 9, 13, 16] },
            r"$B_1^{(1)}$" : {"plus" : [4, 11, 15], "minus" : [1, 8, 12] },
            r"$B_1^{(2)}$" : {"plus" : [3, 5, 9, 10, 16, 17], "minus" : [0, 2, 6, 7, 13, 14] },
            r"$B_2^{(1)}$" : {"plus" : [2, 3, 6, 10, 13, 17], "minus" : [0, 5, 7, 9, 14, 16] },
            r"$E_1^{(1)}$" : {"plus" : [2, 9], "minus" : [5, 6] },
            r"$E_1^{(2)}$" : {"plus" : [0, 10], "minus" : [3, 7] },
            r"$E_1^{(3)}$" : {"plus" : [1, 11], "minus" : [4, 8] },
            r"$E_1^{(4)}$" : {"plus" : [1, 15], "minus" : [4, 12] },
            r"$E_1^{(5)}$" : {"plus" : [2, 16], "minus" : [5, 13] },
            r"$E_1^{(6)}$" : {"plus" : [0, 17], "minus" : [3, 14] },
            r"$E_2^{(1)}$" : {"plus" : [6, 9], "minus" : [2, 5] },
            r"$E_2^{(2)}$" : {"plus" : [7, 10], "minus" : [0, 3] },
            r"$E_2^{(3)}$" : {"plus" : [8, 11], "minus" : [1, 4] },
            r"$E_2^{(4)}$" : {"plus" : [12, 15], "minus" : [1, 4] },
            r"$E_2^{(5)}$" : {"plus" : [13, 16], "minus" : [2, 5] },
            r"$E_2^{(6)}$" : {"plus" : [14, 17], "minus" : [0, 3] },
        },
        "next" : {
            r"$A_1^{(1)}$" : {"plus" : [2, 5, 7, 10, 12, 15], "minus" : [] },
            r"$A_1^{(2)}$" : {"plus" : [0, 1, 3, 4, 6, 8, 9, 11, 13, 14, 16, 17], "minus" : [] },
            r"$A_2^{(1)}$" : {"plus" : [1, 4, 6, 9, 14, 17], "minus" : [0, 3, 8, 11, 13, 16] },
            r"$B_1^{(1)}$" : {"plus" : [5, 7, 15], "minus" : [2, 10, 12] },
            r"$B_1^{(2)}$" : {"plus" : [1, 3, 9, 11, 13, 17], "minus" : [0, 4, 6, 8, 14, 16] },
            r"$B_2^{(1)}$" : {"plus" : [0, 1, 8, 9, 16, 17], "minus" : [3, 4, 6, 11, 13, 14] },
            r"$E_1^{(1)}$" : {"plus" : [4, 9], "minus" : [1, 6] },
            r"$E_1^{(2)}$" : {"plus" : [5, 10], "minus" : [2, 7] },
            r"$E_1^{(3)}$" : {"plus" : [0, 11], "minus" : [3, 8] },
            r"$E_1^{(4)}$" : {"plus" : [2, 15], "minus" : [5, 12] },
            r"$E_1^{(5)}$" : {"plus" : [3, 16], "minus" : [0, 13] },
            r"$E_1^{(6)}$" : {"plus" : [4, 17], "minus" : [1, 14] },
            r"$E_2^{(1)}$" : {"plus" : [6, 9], "minus" : [1, 4] },
            r"$E_2^{(2)}$" : {"plus" : [7, 10], "minus" : [2, 5] },
            r"$E_2^{(3)}$" : {"plus" : [8, 11], "minus" : [0, 3] },
            r"$E_2^{(4)}$" : {"plus" : [12, 15], "minus" : [2, 5] },
            r"$E_2^{(5)}$" : {"plus" : [13, 16], "minus" : [0, 3] },
            r"$E_2^{(6)}$" : {"plus" : [14, 17], "minus" : [1, 4] },
        },
    }

  def getSymmetryNames(self) -> tuple[str, ...]:
    symmetries = (r"A_1^{(1)}", r"A_1^{(2)}",
                  r"A_2^{(1)}",
                  r"B_1^{(1)}", r"B_1^{(2)}",
                  r"B_2^{(1)}",
                  r"E_1^{(1)}", r"E_1^{(2)}", r"E_1^{(3)}", r"E_1^{(4)}", r"E_1^{(5)}", r"E_1^{(6)}",
                  r"E_2^{(1)}", r"E_2^{(2)}", r"E_2^{(3)}", r"E_2^{(4)}", r"E_2^{(5)}", r"E_2^{(6)}")
    return symmetries

  def getProjectionIndeces(self) -> dict[str, dict[str, dict[str,list[int]]]]:
    return self.projectionIndeces

  def getCallbacks(self) -> list[Callable[[list[np.complex128], str], np.complex128]]:
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
  def A1Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$A_1^{(1)}$"]["plus"],
                            self.projectionIndeces[neighbourType][r"$A_1^{(1)}$"]["minus"])

  def A1Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$A_1^{(2)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$A_1^{(2)}$"]["minus"])

  def A2Projection(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$A_2^{(1)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$A_2^{(1)}$"]["minus"])

  # f-wave
  def B1Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$B_1^{(1)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$B_1^{(1)}$"]["minus"])

  def B1Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$B_1^{(2)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$B_1^{(2)}$"]["minus"])


  def B2Projection(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$B_2^{(1)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$B_2^{(1)}$"]["minus"])


  # p-wave
  def E1Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(1)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(1)}$"]["minus"])


  def E1Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(2)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(2)}$"]["minus"])
  def E1Projection3(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(3)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(3)}$"]["minus"])

  def E1Projection4(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(4)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(4)}$"]["minus"])


  def E1Projection5(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(5)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(5)}$"]["minus"])


  def E1Projection6(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_1^{(6)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_1^{(6)}$"]["minus"])

  # d-wave
  def E2Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(1)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(1)}$"]["minus"])

  def E2Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(2)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(2)}$"]["minus"])


  def E2Projection3(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(3)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(3)}$"]["minus"])


  def E2Projection4(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(4)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(4)}$"]["minus"])


  def E2Projection5(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(5)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(5)}$"]["minus"])


  def E2Projection6(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType][r"$E_2^{(6)}$"]["plus"],
                           self.projectionIndeces[neighbourType][r"$E_2^{(6)}$"]["minus"])


class ProjectorC3v(Projector):

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  def __init__(self):
    self.projectionIndeces: dict[str, dict[str, dict[str, list[int]]]] = {
        "nearest" : {
            r"$A_1^{(1)}": {"plus" : [1, 8, 12], "minus": []},
            r"$A_1^{(2)}": {"plus" : [0, 2, 6, 7, 13, 14], "minus": []},
            r"$A_1^{(3)}": {"plus" : [4, 11, 15], "minus": []},
            r"$A_1^{(4)}": {"plus" : [3, 5, 9, 10, 16, 17], "minus": []},
            r"$A_2^{(1)}": {"plus" : [0, 7, 14], "minus": [2, 6, 13]},
            r"$A_2^{(2)}": {"plus" : [3, 10, 17], "minus": [5, 9, 16]},
            r"$E^{(1)}": {"plus" : [6], "minus": [2]},
            r"$E^{(2)}": {"plus" : [7], "minus": [0]},
            r"$E^{(3)}": {"plus" : [8], "minus": [1]},
            r"$E^{(4)}": {"plus" : [9], "minus": [5]},
            r"$E^{(5)}": {"plus" : [10], "minus": [3]},
            r"$E^{(6)}": {"plus" : [11], "minus": [4]},
            r"$E^{(7)}": {"plus" : [12], "minus": [1]},
            r"$E^{(8)}": {"plus" : [13], "minus": [2]},
            r"$E^{(9)}": {"plus" : [14], "minus": [0]},
            r"$E^{(10)}": {"plus" : [15], "minus": [4]},
            r"$E^{(11)}": {"plus" : [16], "minus": [5]},
            r"$E^{(12)}": {"plus" : [17], "minus": [3]},
        },
    }


  def getSymmetryNames(self) -> tuple[str, ...]:
    symmetries = (r"$A_1^{(1)}$", r"$A_1^{(2)}$", r"$A_1^{(3)}$", r"$A_1^{(4)}$",
                  r"$A_2^{(1)}$", r"$A_2^{(2)}$",
                  r"$E^{(1)}$", r"$E^{(2)}$", r"$E^{(3)}$", r"$E^{(4)}$", r"$E^{(5)}$", r"$E^{(6)}$",
                  r"$E^{(7)}$", r"$E^{(8)}$", r"$E^{(9)}$", r"$E^{(10)}$", r"$E^{(11)}$", r"$E^{(12)}$")
    return symmetries

  def getProjectionIndeces(self) -> dict[str, dict[str, dict[str, list[int]]]]:
    return self.projectionIndeces

  def getCallbacks(self) -> list[Callable[[list[np.complex128], str], np.complex128]]:
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
  def A1Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A1_1"]["plus"],
                           self.projectionIndeces[neighbourType]["A1_1"]["minus"])

  def A1Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A1_2"]["plus"],
                           self.projectionIndeces[neighbourType]["A1_2"]["minus"])

  def A1Projection3(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A1_3"]["plus"],
                           self.projectionIndeces[neighbourType]["A1_3"]["minus"])

  def A1Projection4(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A1_4"]["plus"],
                           self.projectionIndeces[neighbourType]["A1_4"]["minus"])

  def A2Projection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A2_1"]["plus"],
                           self.projectionIndeces[neighbourType]["A2_1"]["minus"])

  def A2Projection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["A2_2"]["plus"],
                           self.projectionIndeces[neighbourType]["A2_2"]["minus"])

  def EProjection1(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_1"]["plus"],
                           self.projectionIndeces[neighbourType]["E_1"]["minus"])

  def EProjection2(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_2"]["plus"],
                           self.projectionIndeces[neighbourType]["E_2"]["minus"])


  def EProjection3(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_3"]["plus"],
                           self.projectionIndeces[neighbourType]["E_3"]["minus"])

  def EProjection4(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_4"]["plus"],
                           self.projectionIndeces[neighbourType]["E_4"]["minus"])

  def EProjection5(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_5"]["plus"],
                           self.projectionIndeces[neighbourType]["E_5"]["minus"])

  def EProjection6(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_6"]["plus"],
                           self.projectionIndeces[neighbourType]["E_6"]["minus"])

  def EProjection7(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_7"]["plus"],
                           self.projectionIndeces[neighbourType]["E_7"]["minus"])

  def EProjection8(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_8"]["plus"],
                           self.projectionIndeces[neighbourType]["E_8"]["minus"])

  def EProjection9(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_9"]["plus"],
                           self.projectionIndeces[neighbourType]["E_9"]["minus"])

  def EProjection10(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_10"]["plus"],
                           self.projectionIndeces[neighbourType]["E_10"]["minus"])

  def EProjection11(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_11"]["plus"],
                           self.projectionIndeces[neighbourType]["E_11"]["minus"])

  def EProjection12(self, listOfGammas: list[np.complex128], neighbourType: str) -> np.complex128:
      return self._getProj(listOfGammas,
                           self.projectionIndeces[neighbourType]["E_12"]["plus"],
                           self.projectionIndeces[neighbourType]["E_12"]["minus"])




