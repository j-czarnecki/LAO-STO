from scipy.optimize import differential_evolution
import numpy as np
import os
import csv
import subprocess
import fortranformat as ff
import pandas as pd
from DataReaderClass import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error
import f90nml
import time

class DosFitter():
  def __init__(self):
    self.nEvals = 0
    self.subbands = 1
    self.sublattices = 3
    self.nNeighbors = 3
    self.nNextNeighbors = 6
    self.layerCouplings = 2 * (self.sublattices - 1)
    self.orbitals = 3
    self.dataReader = DataReader("./", "xxx", self.sublattices, self.subbands)
    self.timeStart = 0


  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """
  def fit(self):
    result = differential_evolution(self.__costMse,
                                    bounds = [(0, 0.5), (-0.3, 0.3), (-0.3, 0.3)],
                                    disp = True)
    print(result, flush=True)

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Private methods ------------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  def __costMse(self, x: list, *args) -> float:
    """
    x = [eFermi, gamma1, gamma2, gamma3]
    """

    self.timeStart = time.time()

    dosFit = self.__getDos(x)
    dosExp = pd.read_csv("/home/pwojcik/NicolasFit/J1_exp_0mT.dat", delim_whitespace=True, names=["E", "DOS"], dtype=np.float64)
    dosExp["DOS"] = dosExp["DOS"] / dosExp["DOS"].max() #Normalization to maximum
    dosExp["E"] = dosExp["E"] * 1000 #Conversion to meV
    dosExp.drop(dosExp[np.abs(dosExp["E"]) > 0.35].index, inplace=True)

    interpolate = interp1d(dosFit["E"], dosFit["DOS"], kind="linear", fill_value="extrapolate")
    DosInterpolated = interpolate(dosExp["E"]) #Interpolate fited DOS at those points where experimental measurements were made

    weightsExponent = 3.0
    error = mean_squared_error(dosExp["DOS"], DosInterpolated)
    print("MSE: ", error, flush=True)

    plt.figure()
    plt.plot(dosFit["E"], dosFit["DOS"], label="fit")
    plt.plot(dosExp["E"], dosExp["DOS"], label="exp")
    plt.plot(dosExp["E"], DosInterpolated, label="interpolated", linestyle="dashed")
    plt.legend()
    plt.savefig(f"../Plots/DosFitHistory_{self.nEvals % 5}.png")
    plt.close()

    self.nEvals += 1
    print(f"x = {x}", flush=True)
    print(f"Evaluation {self.nEvals} took: {time.time() - self.timeStart} seconds.", flush=True)

    return error

  def __getDos(self, x: list) -> pd.DataFrame:
    eFermi = x[0]
    gammas = x[1:]
    gammas = [gammas[0], gammas[0], gammas[1]] #Assume that coupling in two directions is the same and only one is different - as in s-wave from self-consistency
    self.__constructAndWriteGamma(gammas)
    self.__changeFermiEnergy(eFermi)
    os.chdir("/home/pwojcik/LAO-STO")
    result = subprocess.run("srun -c 64 bin/POST_LAO_STO.x", shell=True, check=True)
    os.chdir("/home/pwojcik/LAO-STO/Analyzer")
    self.dataReader.LoadDos("/home/pwojcik/LAO-STO/OutputData/DOS.dat")
    result = subprocess.run(f"cp /home/pwojcik/LAO-STO/OutputData/DOS.dat /home/pwojcik/DOS_train/DOS_{self.nEvals}.dat", shell=True, check=True)
    with open(f"/home/pwojcik/DOS_train/X_{self.nEvals}.dat", "w") as file:
      print(f"{x[0]} {x[1]} {x[2]}", file=file)
    return self.dataReader.dosDataframe

  def __constructAndWriteGamma(self, gammas: list) -> None:
    if len(gammas) != 3:
      raise ValueError("List of gammas must have 3 elements")

    neigh_orb_start_map = {1: 0, 2: 2, 3: 1}

    with open("../OutputData/Gamma_SC_final.dat", "w", newline="") as file:
      print(" #band spin neighbour lattice orbital Re(Gamma) Im(Gamma)", file=file)
      fortFormat = ff.FortranRecordWriter("(5I5, 2E15.5)")
      for band in range(1, max(1, self.subbands) + 1):
        for spin in range(1, 3):
          spinSign = (-1)**(spin - 1)
          for neigh in range(1, self.nNeighbors + self.nNextNeighbors + 1):
              maxLat = self.layerCouplings + 1 if neigh <= self.nNeighbors else self.sublattices + 1
              for lat in range(1, maxLat):
                for orb in range(1, self.orbitals + 1):
                  if neigh <= self.nNeighbors:
                    gamma_idx = (neigh_orb_start_map[neigh] + orb - 1) % self.orbitals
                    line = fortFormat.write([band, spin, neigh, lat, orb, spinSign * gammas[gamma_idx], 0])
                    print(line, file=file)
                  else:
                    gamma_idx = (neigh - 1 + orb - 1) % self.orbitals
                    line = fortFormat.write([band, spin, neigh, lat, orb, 0, 0])
                    print(line, file=file)

              print(" ", file=file)
              print(" ", file=file)

  def __changeFermiEnergy(self, energy: float) -> None:
    nml = f90nml.read("/home/pwojcik/LAO-STO/input.nml")
    nml["physical_params"]["E_Fermi"] = energy * 1000 # Convert eV to meV
    nml.write("/home/pwojcik/LAO-STO/input.nml", force = True)

def main():
  dosFit = DosFitter()
  dosFit.fit()



if __name__ == "__main__":
  main()








# def costFunction(x: list, *args):
#   return np.exp(-x[0]**2)*(x[0] - 2)**3 + (x[1] + 3)**2 + (x[2] - 1)**2 + np.exp(x[0])



# bounds = [(-5, 5) for i in range(3)]
# print(bounds)
# for i in range(5):
#   result = differential_evolution(costFunction, bounds)
#   print(result.x[0])
#   print(f"x = {1 - np.sqrt(5/2)} y = {-3} z = {1}")

#   print(costFunction([1 - np.sqrt(5/2), -3, 1], 0))