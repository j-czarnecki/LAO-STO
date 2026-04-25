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

import logging
import sys

import sympy as sp

# caution: path[0] is reserved for script path (or '' in REPL)
from GroupTheory.GroupGeneratorsClass import C3vGenerator
from GroupTheory.SymbolicSymmetryProjectorClass import SymbolicSymmetryProjectorClass
from sympy.physics.quantum import TensorProduct

logger = logging.getLogger(__name__)


class SymmetryResolver:
  def __init__(self, nSublats: int):
    """
    Parameters
    ----------
    nSublats : int
      Number of sublattices in the system
    """
    self.nSublats: int = nSublats

    # Real-space submodule
    # On-site pairing
    self.c3Onsite: sp.Matrix
    self.sV1Onsite: sp.Matrix
    # Nearest-neighbours pairing
    self.c3Nearest: sp.Matrix
    self.sV1Nearest: sp.Matrix
    # Next-neighbours pairing
    self.c3Next: sp.Matrix
    self.sV1Next: sp.Matrix

    # Orbital submodule
    self.c3Orb: sp.Matrix
    self.sV1Orb: sp.Matrix
    # Spin submodule
    self.c3Spin: sp.Matrix
    self.sV1Spin: sp.Matrix

    # Basis pairing matrices
    # Orbital submodule
    # Intraorbital pairing
    self.orbYz: sp.Matrix
    self.orbZx: sp.Matrix
    self.orbXy: sp.Matrix
    # Even interorbital pairing
    self.orbEvenYzZx: sp.Matrix
    self.orbEvenYzXy: sp.Matrix
    self.orbEvenZxXy: sp.Matrix
    # Odd interorbital pairing
    self.orbOddYzZx: sp.Matrix
    self.orbOddYzXy: sp.Matrix
    self.orbOddZxXy: sp.Matrix
    self.orbBasis: list[sp.Matrix]

    # Sublattice submodule
    # TODO: think about the construction that will take into account multi-lattice (>2) setups
    self.interLat: sp.Matrix
    self.intraLat: sp.Matrix
    self.interLatBasis: sp.Matrix
    self.intraLatBasis: sp.Matrix

    # Spin submodule
    # Pauli matrices
    self.s0: sp.Matrix
    self.sx: sp.Matrix
    self.sy: sp.Matrix
    self.sz: sp.Matrix
    self.spinBasis: list[sp.Matrix]

    self.kx, self.ky = sp.symbols("kx ky", real=True)
    self.k1 = self.ky
    self.k2 = sp.sqrt(3) / sp.S(2) * self.kx - sp.Rational(1, 2) * self.ky
    self.k3 = -sp.sqrt(3) / sp.S(2) * self.kx - sp.Rational(1, 2) * self.ky

    self.k1Next = sp.sqrt(3) * self.kx
    self.k2Next = sp.sqrt(3) / sp.S(2) * self.kx + sp.Rational(3, 2) * self.ky
    self.k3Next = -sp.sqrt(3) / sp.S(2) * self.kx + sp.Rational(3, 2) * self.ky

    self.__createGenerators()
    self.__createBasisPairingMatrices()

  def __createGenerators(self):
    """
    Initializes generators of symmetry in all submodules.
    """
    # Real-space submodule
    # On-site pairing
    self.c3Onsite = sp.Matrix([[1]])
    self.sV1Onsite = sp.Matrix([[1]])
    # Nearest-neighbours pairing
    self.c3Nearest = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    self.sV1Nearest = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    # Next-neighbors pairing
    self.c3Next = sp.Matrix([
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1],
      [1, 0, 0, 0, 0, 0],
    ])
    self.sV1Next = sp.Matrix([
      [1, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
    ])

    # Orbital submodule
    self.c3Orb = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    self.sV1Orb = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    # Spin submodule
    self.c3Spin = sp.Matrix([
      [1, 0, 0, 0],
      [0, sp.cos(sp.pi * 2 / 3), sp.sin(-sp.pi * 2 / 3), 0],
      [0, sp.sin(sp.pi * 2 / 3), sp.cos(sp.pi * 2 / 3), 0],
      [0, 0, 0, 1],
    ])
    self.sV1Spin = sp.Matrix([
      [1, 0, 0, 0],
      [0, sp.Rational(1, 2), sp.sqrt(3) / sp.S(2), 0],
      [0, sp.sqrt(3) / sp.S(2), -sp.Rational(1, 2), 0],
      [0, 0, 0, 1],
    ])

  def __createBasisPairingMatrices(self):
    """
    Construct basis of pairing matrices for each submodule.
    """
    # Define matrices that represent intraorbital coupling for a given index
    # in the irreducible representation eigenvector.
    # Intraorbital coupling
    self.orbYz = sp.Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
    self.orbZx = sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    self.orbXy = sp.Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]])

    # Interorbital-even coupling
    self.orbEvenYzZx = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    self.orbEvenYzXy = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    self.orbEvenZxXy = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])

    # Interorbital-odd coupling
    self.orbOddYzZx = sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
    self.orbOddYzXy = sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]])
    self.orbOddZxXy = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
    self.orbBasis = [
      self.orbYz,
      self.orbZx,
      self.orbXy,
      self.orbEvenYzXy,
      self.orbEvenYzXy,
      self.orbEvenZxXy,
      self.orbOddYzXy,
      self.orbOddYzXy,
      self.orbOddZxXy,
    ]

    self.interLat = sp.zeros(self.nSublats)
    for i in range(self.nSublats - 1):
      self.interLat[i, i + 1] = 1  # Right stripe (upper diagonal)
      self.interLat[i + 1, i] = 1  # Left stripe (lower diagonal)

    self.interLatBasis = sp.Matrix([sp.exp(sp.I * self.k3), sp.exp(sp.I * self.k2), sp.exp(sp.I * self.k1)])
    self.intraLatBasis = sp.Matrix([1, self.k1Next, self.k2Next, self.k3Next, -self.k1Next, -self.k2Next, -self.k3Next])

    # Pauli matrices
    self.s0 = sp.eye(2)
    self.sx = sp.Matrix([[0, 1], [1, 0]])
    self.sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    self.sz = sp.Matrix([[1, 0], [0, -1]])
    self.spinBasis = [self.s0, self.sx, self.sy, self.sz]

  def __getOrbChannelToCanonicalTransformation(self) -> sp.Matrix:
    """
    Get transformation between canonical ordering of orbital pairing matrix
    and channel-resolved ordering. In the latter, the basis matrices for orbital pairing are
    intraorbital, even interorbital and odd interorbital.
    """
    # Create 9x9 matrix, where each column corresponds to vectorized 3x3 basis matrix M
    nOrbElems = self.orbYz.shape[0] * self.orbYz.shape[1]
    B = sp.Matrix.hstack(*[M.reshape(nOrbElems, 1) for M in self.orbBasis])
    return B

  def __getProjectionsEigenvectors(self) -> dict[str, list[tuple[int, int, list[sp.Matrix]]]]:
    """
    Solve eigenproblem for projection operators and return set of eigenvectors corresponding
    to eigenvalues 1 (projected states) for each irreducible representation.
    """
    # We have to take into consideration all possible hoppings
    c3NeighborsTotal = sp.Matrix.diag(self.c3Onsite, self.c3Nearest, self.c3Nearest)
    sV1NeighborsTotal = sp.Matrix.diag(self.sV1Onsite, self.sV1Nearest, self.sV1Nearest)

    # Take tensor product of orbital representations to account for
    # intra- and interorbital pairing.
    # Apply basis change B, since we want to move from canonical ordering (d_11, d_12, d_13, d_21, ..., d_33)
    # To a new basis of pairing channels.
    B = self.__getOrbChannelToCanonicalTransformation()
    c3OrbTotal = B.inv() * TensorProduct(self.c3Orb, self.c3Orb) * B
    sV1OrbTotal = B.inv() * TensorProduct(self.sV1Orb, self.sV1Orb) * B

    c3Total = TensorProduct(self.c3Spin, c3OrbTotal, c3NeighborsTotal)
    sV1Total = TensorProduct(self.sV1Spin, sV1OrbTotal, sV1NeighborsTotal)

    c3vGroup = C3vGenerator(c3Total, sV1Total)

    symmetryResolver = SymbolicSymmetryProjectorClass(
      c3vGroup.irrepsTuple, c3vGroup.conjugacyClassesTuple, c3vGroup.chiTabDict, c3vGroup.representationDim
    )
    projectionMatrices = symmetryResolver.getProjectionOperators(c3vGroup.getOperationsDict())
    multiplicities = symmetryResolver.getMultiplicities(c3vGroup.getOperationsDict())
    for irrep in multiplicities:
      logger.info(f"N({irrep}) = {multiplicities[irrep]}")
    eigenproblem = symmetryResolver.getDiagonalizedProjections(projectionMatrices)
    eigenproblem = symmetryResolver.filterOutProjections(eigenproblem)
    return eigenproblem

  def __getSpinAndOrbitalIndicesFromVectorizedMatrix(self, index) -> tuple[int, int]:
    """
    Get spin and orbital indices of a given vectorized matrix.
    The assumption is that the matrix has already be reshaped so that each element of a vector
    contains a vector with neighbors pairing and the only indices left are spin and orbital.

    Parameters:
    index: int
        Index in spin-orbital vectorized matrix

    Returns:
    spinIndex, orbIndex: int, int
        Index of spin and orbital basis matrix
    """
    nOrbBasisMatrices = len(self.orbBasis)

    spin = index // nOrbBasisMatrices
    orb = index % nOrbBasisMatrices

    return spin, orb

  def __matricizeIrrepEigenvector(self, eigenvector: sp.Matrix) -> sp.Matrix:
    nOrbBasisMatrices = len(self.orbBasis)
    nSpinBasisMatrices = len(self.spinBasis)
    nNearestNeighborsPerSite = 3
    nNextNeighborsPerSite = 6
    eigenvectorSpinOrbital = eigenvector.reshape(
      nOrbBasisMatrices * nSpinBasisMatrices, nNearestNeighborsPerSite + nNextNeighborsPerSite + 1
    )  # +1 due to onsite pairing

    # nSpin * nOrbitals * nSublattices is the dimension of pairing matrix
    irrepMatrix = sp.zeros((2 * 3 * self.nSublats, 2 * 3 * self.nSublats))

    for i in range(len(eigenvectorSpinOrbital[:, 0])):
      spin, orb = self.__getSpinAndOrbitalIndicesFromVectorizedMatrix(i)
      # Nearest-neighbours pairing
      interLatNeighbors = eigenvectorSpinOrbital[i, 1 : (1 + nNearestNeighborsPerSite)]
      # On-site pairing + next-neighbours pairing
      intraLatNeighbors = sp.Matrix([
        eigenvectorSpinOrbital[i, 0],
        *eigenvectorSpinOrbital[i, (1 + nNearestNeighborsPerSite) :],
      ])

      kElemInterLat = interLatNeighbors.dot(self.interLatBasis)
      kElemIntraLat = intraLatNeighbors.dot(self.intraLatBasis)

      elemMatrixInterLat = TensorProduct(-sp.I * self.spinBasis[spin] * self.sy, self.interLat, self.orbBasis[orb])
      elemMatrixIntraLat = TensorProduct(-sp.I * self.spinBasis[spin] * self.sy, self.intraLat, self.orbBasis[orb])

      elemMatrixInterLat *= kElemInterLat
      elemMatrixIntraLat *= kElemIntraLat

      irrepMatrix += elemMatrixInterLat
      irrepMatrix += elemMatrixIntraLat

    # Impose fermionic antisymmetry on sublattice hoppings
    irrepMatrixMinusK = irrepMatrix.subs([(self.kx, -self.kx), (self.ky, -self.ky)])

    irrepMatrixAntisymmetrized = irrepMatrix - irrepMatrixMinusK

    return irrepMatrixAntisymmetrized
