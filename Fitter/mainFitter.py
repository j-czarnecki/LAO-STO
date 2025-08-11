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

from Fitter.DosFitter import *
import argparse
import yaml
import logging

def loadConfig(configPath):
  with open(configPath, 'r') as f:
    return yaml.safe_load(f)

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", type=str, default="Fitter/fitConfig.yaml", help="Path to config file")
  args = parser.parse_args()

  print("Args parsed", flush=True)

  config = loadConfig(args.config)
  print("Config read", config, flush=True)


  dosFit = DosFitter(runsDir=config["runsDir"],
                     dosExpPath=config["dosExpPath"],
                     eMax=config["eMax"])
  dosFit.fit(paramBounds=config["bounds"])

if __name__ == "__main__":
  main()
