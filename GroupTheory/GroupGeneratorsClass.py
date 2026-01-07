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


class C6vGenerator:

  def __init__(self, c6: sp.Matrix, s_v1: sp.Matrix, s_d1: sp.Matrix):
    """
    Initializes C6vGenerator object to hold all C6v point group data: irreps, conjugacy classes, character table and generators.
    """
    self.representationDim = c6.shape[0]
    self.c6: sp.Matrix = c6
    self.s_v1: sp.Matrix = s_v1
    self.s_d1: sp.Matrix = s_d1
    self.irrepsTuple = (r"$A_1$", r"$A_2$", r"$B_1$", r"$B_2$", r"$E_1$", r"$E_2$")
    self.conjugacyClassesTuple = ("E", "2C_6", "2C_3", "C_2", "3s_v", "3s_d")
    self.chiTabDict = {r"$A_1$": {"E" : sp.S(1), "2C_6" : sp.S(1) , "2C_3" : sp.S(1) , "C_2" : sp.S(1) , "3s_v" : sp.S(1) , "3s_d" : sp.S(1)},
                       r"$A_2$": {"E" : sp.S(1), "2C_6" : sp.S(1) , "2C_3" : sp.S(1) , "C_2" : sp.S(1) , "3s_v" : sp.S(-1), "3s_d" : sp.S(-1)},
                       r"$B_1$": {"E" : sp.S(1), "2C_6" : sp.S(-1), "2C_3" : sp.S(1) , "C_2" : sp.S(-1), "3s_v" : sp.S(1) , "3s_d" : sp.S(-1)},
                       r"$B_2$": {"E" : sp.S(1), "2C_6" : sp.S(-1), "2C_3" : sp.S(1) , "C_2" : sp.S(-1), "3s_v" : sp.S(-1), "3s_d" : sp.S(1)},
                       r"$E_1$": {"E" : sp.S(2), "2C_6" : sp.S(1) , "2C_3" : sp.S(-1), "C_2" : sp.S(-2), "3s_v" : sp.S(0) , "3s_d" : sp.S(0)},
                       r"$E_2$": {"E" : sp.S(2), "2C_6" : sp.S(-1), "2C_3" : sp.S(-1), "C_2" : sp.S(2) , "3s_v" : sp.S(0) , "3s_d" : sp.S(0)}}

  def getOperationsDict(self) -> dict[str, list[sp.Matrix]]:
    e = sp.eye(self.representationDim)
    c6_inv = self.c6.inv()

    c3 = self.c6**2
    c3_inv = c3.inv()

    c2 = self.c6**3

    s_v2 = c6_inv * self.s_v1 * self.c6
    s_v3 = c3_inv * self.s_v1 * c3

    s_d2 = c6_inv * self.s_d1 * self.c6
    s_d3 = c3_inv * self.s_d1 * c3
    return {"E": [e], "2C_6": [self.c6, c6_inv], "2C_3": [c3, c3_inv], "C_2": [c2], "3s_v": [self.s_v1, s_v2, s_v3], "3s_d": [self.s_d1, s_d2, s_d3]}

class C3vGenerator:

  def __init__(self, c3: sp.Matrix, s_v1: sp.Matrix):
    """
    Initializes C6vGenerator object to hold all C6v point group data: irreps, conjugacy classes, character table and generators.
    """
    self.representationDim = c3.shape[0]
    self.c3: sp.Matrix = c3
    self.s_v1: sp.Matrix = s_v1
    self.irrepsTuple = ("A_1", "A_2", "E")
    self.conjugacyClassesTuple = ("E", "2C_3", "3s_v")
    self.chiTabDict = {"A_1": {"E" : sp.S(1), "2C_3" : sp.S(1) , "3s_v" : sp.S(1)},
                       "A_2": {"E" : sp.S(1), "2C_3" : sp.S(1) , "3s_v" : sp.S(-1)},
                       "E":   {"E" : sp.S(2), "2C_3" : sp.S(-1) , "3s_v" : sp.S(0)}}

  def getOperationsDict(self) -> dict[str, list[sp.Matrix]]:
    e = sp.eye(self.representationDim)

    c3_inv = self.c3.inv()

    s_v2 = c3_inv * self.s_v1 * self.c3
    s_v3 = c3_inv * s_v2 * self.c3

    return {"E": [e], "2C_3": [self.c3, c3_inv], "3s_v": [self.s_v1, s_v2, s_v3]}