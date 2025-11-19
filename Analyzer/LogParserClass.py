import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

class LogParser:
  def __init__(self):
    self.iterRegex = re.compile(r"SC_ITER:\s*(\d+)")
    self.divergentIntegrationRegex = re.compile(
        r"TID\s*=\s*(\d+).*?"
        r"k1_chunk_min:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)"
        r"\s*k2_chunk_min:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)"
        r".*?max_error_delta:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)"
        r"\s*max_error_charge:\s*([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)"
    )

  def getDivergentChunks(self, path):
    results = []
    currentIter = None

    with open(path, "r") as f:
      for line in f:
        iterationMatch = self.iterRegex.search(line)
        if iterationMatch:
          currentIter = int(iterationMatch.group(1))
          continue

        divergentMatch = self.divergentIntegrationRegex.search(line)
        if divergentMatch and currentIter is not None:
          tid = int(divergentMatch.group(1))
          k1 = float(divergentMatch.group(2))
          k2 = float(divergentMatch.group(3))
          maxDeltaError = float(divergentMatch.group(4))
          maxChargeError = float(divergentMatch.group(5))
          results.append({"iter": currentIter,
                          "tid": tid,
                          "k1": k1,
                          "k2": k2,
                          "max_delta_error": maxDeltaError,
                          "max_charge_error": maxChargeError})

    df = pd.DataFrame(results)
    df["kx"] = df["k1"] * np.cos(df["k2"])
    df["ky"] = df["k1"] * np.sin(df["k2"])
    return df

  def plotDivergentChunks(self, chunksDf):
    groups = chunksDf.groupby("iter")
    errorTypes = ["delta", "charge"]

    for i, group in groups:
      fig = plt.figure(figsize=(7, 5), dpi=400)
      gs = gridspec.GridSpec(1, 1, figure=fig, left=0.25, right=0.9, top=0.95, bottom=0.1)
      ax = fig.add_subplot(gs[0,0])

      self.__plotFirstBrillouinZoneBoundary()
      scatter = ax.scatter(group["kx"], group["ky"], c=group["max_delta_error"], norm=LogNorm(), marker=".", s=1)
      cbar = plt.colorbar(scatter, ax=ax)
      cbar.set_label("Max delta error")

      ax.set_title("Divergent chunks")
      ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
      ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")

      plt.grid(True)
      plt.savefig(f"../Plots/Divergence_{i}.png")

  def __plotFirstBrillouinZoneBoundary(self, ax = None):
      brillouinZoneVertices = np.zeros((7, 2))  # One more to close the polygon

      brillouinZoneVertices[:, 0] = np.array(
          [
              4.0 * np.pi / (3 * np.sqrt(3.0)),
              2.0 * np.pi / (3 * np.sqrt(3.0)),
              -2.0 * np.pi / (3 * np.sqrt(3.0)),
              -4.0 * np.pi / (3 * np.sqrt(3.0)),
              -2.0 * np.pi / (3 * np.sqrt(3.0)),
              2.0 * np.pi / (3 * np.sqrt(3.0)),
              4.0 * np.pi / (3 * np.sqrt(3.0)),
          ]
      )

      brillouinZoneVertices[:, 1] = np.array(
          [
              0.0,
              -2.0 * np.pi / 3.0,
              -2.0 * np.pi / 3.0,
              0.0,
              2.0 * np.pi / 3.0,
              2.0 * np.pi / 3.0,
              0.0,
          ]
      )
      if ax == None:
        ax = plt.gca()
      ax.plot(
          brillouinZoneVertices[:, 0],
          brillouinZoneVertices[:, 1],
          "--",
          color="black",
          linewidth=2,
      )