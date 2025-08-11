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

from scipy.optimize import differential_evolution
import numpy as np
import os
import csv
import subprocess
import fortranformat as ff
import pandas as pd
from Analyzer.DataReaderClass import *
from Runner.RunnerConfigClass import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error
import f90nml
import time
import seaborn as sns
import argparse
import yaml
import logging

SCRATCH_PATH = os.getenv('SCRATCH')
HOME_PATH = os.getenv('HOME')

class DosFitter():
  def __init__(self, runsDir: str, dosExpPath: str, eMax: float):
    self.nEvals = 0
    self.subbands = 1
    self.sublattices = 3
    self.nNeighbors = 3
    self.nNextNeighbors = 6
    self.layerCouplings = 2 * (self.sublattices - 1)
    self.orbitals = 3
    self.dataReader = DataReader("./", "xxx", self.sublattices, self.subbands)
    self.runnerConfig = RunnerConfig()
    self.timeStart = 0
    self.__initializePlotParams()

    self.runsDir = runsDir
    self.dosExpPath = dosExpPath
    self.eMax = eMax
    self.dosExp = None

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Interface methods ----------------------------------- """
  """ ---------------------------------------------------------------------------------- """
  def fit(self, paramBounds: list[tuple[float, float]]) -> None:

    self.__getDosExp()

    result = differential_evolution(self.__costMse,
                                    bounds = paramBounds,
                                    disp = True)
    print(result, flush=True)

  """ ---------------------------------------------------------------------------------- """
  """ ---------------------------- Private methods ------------------------------------- """
  """ ---------------------------------------------------------------------------------- """

  def __costMse(self, x: list, *args) -> float:
    """
    x = [gamma1, gamma2, gamma3]
    """

    self.timeStart = time.time()

    dosFit = self.__getDos(x)

    # Get appropriate range of DOS and renormalize
    dosFit.drop(dosFit[np.abs(dosFit["E"]) > self.eMax].index, inplace=True)
    dosMax = dosFit["DOS"].max()
    if dosMax == 0.0:
      return np.inf
    dosFit["DOS"] = dosFit["DOS"] / dosMax #Normalization to maximum

    interpolate = interp1d(dosFit["E"], dosFit["DOS"], kind="linear", fill_value="extrapolate")
    DosInterpolated = interpolate(self.dosExp["E"]) #Interpolate fited DOS at those points where experimental measurements were made

    error = mean_squared_error(self.dosExp["DOS"], DosInterpolated)
    print("MSE: ", error, flush=True)

    plt.figure()
    plt.plot(dosFit["E"], dosFit["DOS"], label="Theory", color="black")
    plt.scatter(self.dosExp["E"], self.dosExp["DOS"], label="Exp.", color="red", marker=".")
    plt.xlabel(r"$E~(meV)$")
    plt.ylabel(r"$G~(a.u.)$")
    plt.legend()
    plt.savefig(f"{os.path.join(SCRATCH_PATH, self.runsDir)}/Plots/DosFitHistory_{self.nEvals % 5}.png")
    plt.close()

    print(f"x = {x}", flush=True)
    print(f"Evaluation {self.nEvals} took: {time.time() - self.timeStart} seconds.", flush=True)
    self.nEvals += 1

    return error

  def __getDos(self, x: list) -> pd.DataFrame:
    #eFermi = x[0]
    gammas = x[:] #Got rid of Fermi energy
    gammas = [gammas[0], gammas[1], gammas[2]] #Assume that coupling in two directions is the same and only one is different - as in s-wave from self-consistency
    self.__constructAndWriteGamma(gammas)
    #self.__changeFermiEnergy(eFermi)
    os.chdir(os.path.join(SCRATCH_PATH, self.runsDir))
    result = subprocess.run(f"{os.path.join(HOME_PATH,'LAO-STO','bin','POST_LAO_STO.x')}", shell=True, check=True)
    self.dataReader.LoadDos(os.path.join(SCRATCH_PATH,self.runsDir, "OutputData","DOS.dat"))
    result = subprocess.run(f"cp {os.path.join(SCRATCH_PATH, self.runsDir, 'OutputData', 'DOS.dat')} {os.path.join(SCRATCH_PATH, self.runsDir, 'DOS_train', f'DOS_{self.nEvals}.dat')}", shell=True, check=True)
    with open(f"{os.path.join(SCRATCH_PATH, self.runsDir, 'DOS_train', f'X_{self.nEvals}.dat')}", "w") as file:
      print(f"{x}", file=file)
    return self.dataReader.dosDataframe

  def __getDosExp(self) -> None:
    #Load experimental DOS only once
    self.dosExp = pd.read_csv(self.dosExpPath, delim_whitespace=True, names=["E", "DOS"], dtype=np.float64)
    self.dosExp["E"] = self.dosExp["E"] * 1000 #Conversion to meV
    self.dosExp.drop(self.dosExp[np.abs(self.dosExp["E"]) > self.eMax].index, inplace=True)
    self.dosExp["DOS"] = self.dosExp["DOS"] - self.dosExp["DOS"].min() # Shift so that minimum is at zero
    self.dosExp["DOS"] = self.dosExp["DOS"] / self.dosExp["DOS"].max() #Normalization to maximum

  def __constructAndWriteGamma(self, gammas: list) -> None:
    if len(gammas) != 3:
      raise ValueError("List of gammas must have 3 elements")

    neigh_orb_start_map = {1: 0, 2: 2, 3: 1}

    with open(f"{os.path.join(SCRATCH_PATH, self.runsDir, 'OutputData', 'Gamma_SC_final.dat')}", "w", newline="") as file:
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
    nml = f90nml.read("/home/jczarnecki/LAO-STO/input.nml")
    nml["physical_params"]["E_Fermi"] = energy * 1000 # Convert eV to meV
    nml.write("/home/jczarnecki/LAO-STO/input.nml", force = True)

  def __initializePlotParams(self):
      plt.rcParams["text.usetex"] = True
      plt.rcParams["font.family"] = "serif"
      plt.rcParams["font.serif"] = "Computer Modern Roman"
      plt.rcParams["font.sans-serif"] = "Computer Modern Sans serif"
      plt.rcParams["font.monospace"] = "Computer Modern Typewriter"
      plt.rcParams["axes.titlesize"] = 30
      plt.rcParams["axes.labelsize"] = 30
      plt.rcParams["xtick.labelsize"] = 26
      plt.rcParams["ytick.labelsize"] = 26
      plt.rcParams["legend.fontsize"] = 20
      plt.rcParams["legend.title_fontsize"] = 24
      # Optionally, add custom LaTeX preamble
      plt.rcParams["text.latex.preamble"] = (
          r"\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}"
      )

      self.__setPalette()

      # Set rcParams for tighter layout
      plt.rcParams["figure.autolayout"] = True
      plt.rcParams["figure.constrained_layout.use"] = False
      plt.rcParams["axes.linewidth"] = 1.2

      # Set rcParams to show ticks on both left and right sides
      plt.rcParams["xtick.direction"] = "in"
      plt.rcParams["ytick.direction"] = "in"
      plt.rcParams["xtick.bottom"] = True
      plt.rcParams["ytick.left"] = True
      plt.rcParams["xtick.top"] = True
      plt.rcParams["ytick.right"] = True

      plt.rcParams["axes.xmargin"] = 0.01

  def __setPalette(self, nColors: int = 3, palette: str = "colorblind"):
      # Choose a seaborn palette
      # has to specify number of lines
      self.palette = sns.color_palette(palette, nColors)
      # Set the color cycle
      plt.rcParams["axes.prop_cycle"] = plt.cycler(color=self.palette)


