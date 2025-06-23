import sympy as sp
import numpy as np
from sympy.vector import CoordSys3D
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import PowerNorm, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from IPython.display import display
import sys
import os
from contextlib import contextmanager

try:
    from IPython.utils.io import capture_output
except ImportError:
    capture_output = None  # fallback if not in IPython


@contextmanager
def initPrinting(enablePrinting=True):
  """ Initializes sympy printing, or disables printing at all"""
  if not enablePrinting:
    savedStdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
       if capture_output is not None:
         with capture_output():
           yield
       else:
         yield
    finally:
      sys.stdout.close()
      sys.stdout = savedStdout
  else:
    sp.init_printing()
    yield


class SymmetryPlotterClass:
  def __init__(self, gridPoints: int = 1000):
    self.__setPlotParams()
    self.nCalls = 0
    self.gridPoints = gridPoints
    self.kMesh = self.__createKMesh()
    # self.cartesianKSymbols
    # self.grapheneKSymbols
    # self.orbitalShapeSymbols
    # self.substitutionsNearest
    # self.substitutionsNext
    self.__createSymbols()
    self.__createSubstitutions()
    #self.fOrbitals
    self.__createLamdifiedOrbitalShapeFunctions()
    self.__createCoordinateSystemAndDeltaTabs()

  def plotBasisFunctions(self, irrepsEigenvectors: dict[str, list[tuple[int, int, list[sp.Matrix]]]], enablePrinting=True):
    with initPrinting(enablePrinting):
      for irrep in irrepsEigenvectors:
        print(irrep)

        for v in irrepsEigenvectors[irrep]:
          if np.abs(v[0]) == 1:
            for i in range(v[1]):
              self.__getSingleBasisFunction(v[2][i])

        print("------------------------------------------------------------")
        print("\n")

  def plotOrbitalWeights(self):
    f_yz, f_zx, f_xy = self.fOrbitals
    X, Y = self.kMesh
    #self.__plotRGBLegend()

    #Orbital shape functions
    R_yz = f_yz(X[:-1, :-1], Y[:-1, :-1])
    G_zx = f_zx(X[:-1, :-1], Y[:-1, :-1])
    B_xy = f_xy(X[:-1, :-1], Y[:-1, :-1])
    #Shift to minimal values equal 0
    R_yz = R_yz - np.min(R_yz)
    G_zx = G_zx - np.min(G_zx)
    B_xy = B_xy - np.min(B_xy)
    #Normalize
    R_yz = R_yz / np.max(R_yz)
    G_zx = G_zx / np.max(G_zx)
    B_xy = B_xy / np.max(B_xy)
    RGB = np.stack([R_yz, G_zx, B_xy], axis=-1).reshape(-1,3)

    fig = plt.figure(figsize=(7, 5), dpi=100)
    # Set up GridSpec (1 row, 1 column, with some spacing)
    gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
    ax = fig.add_subplot(gs[0,0])
    cax = ax.pcolormesh(X, Y, np.zeros_like(R_yz), color=RGB, shading='flat')
    self.__plotFirstBrillouinZoneBoundary()
    ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
    ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.show()
    plt.close()


  def __createKMesh(self):
    kx = np.linspace(-2.5, 2.5, self.gridPoints)
    ky = np.linspace(-2.5, 2.5, self.gridPoints)
    return np.meshgrid(kx, ky)

  def __createSymbols(self):
    self.cartesianKSymbols = sp.symbols(r'k_x, k_y, k_z')
    self.grapheneKSymbols = sp.symbols(r'k_1, k_2, k_3')
    self.orbitalShapeSymbols = sp.symbols(r'c_yz, c_zx, c_xy')

  def __createSubstitutions(self):
    k_x, k_y, _ = self.cartesianKSymbols
    k_1, k_2, k_3 = self.grapheneKSymbols
    self.substitutionsKGraphene = ((k_1, sp.sqrt(3)/sp.S(2) * k_x + sp.Rational(1,2) * k_y),
                                   (k_2, -sp.sqrt(3)/sp.S(2) * k_x + sp.Rational(1,2) * k_y),
                                   (k_3, -k_y))
    self.substitutionsNearest = ((sp.sqrt(3)/sp.S(2) * k_x + sp.Rational(1,2) * k_y, k_1),
                                 (-sp.sqrt(3)/sp.S(2) * k_x + sp.Rational(1,2) * k_y, k_2),
                                 (k_y, -k_3),
                                 (sp.sqrt(3)/sp.S(2) * k_x, sp.Rational(1,2) * (k_1 - k_2)))
    self.substitutionsNext = ((-sp.sqrt(3)/sp.S(2) * k_x + sp.Rational(3,2) * k_y, k_1),
                              (-sp.sqrt(3)/sp.S(2) * k_x - sp.Rational(3,2) * k_y, k_2),
                              (sp.sqrt(3) * k_x, k_3),
                              (sp.Rational(3,2) * k_y, sp.Rational(1,2) * (k_1 - k_2)))

  def __createLamdifiedOrbitalShapeFunctions(self):
    k_x, k_y, _ = self.cartesianKSymbols
    c_yz, c_zx, c_xy = self.__getSymbolicOrbitalShapeFunctions()
    f_yz = sp.lambdify((k_x, k_y), c_yz, "numpy")
    f_zx = sp.lambdify((k_x, k_y), c_zx, "numpy")
    f_xy = sp.lambdify((k_x, k_y), c_xy, "numpy")
    self.fOrbitals = (f_yz, f_zx, f_xy)

  def __getSymbolicOrbitalShapeFunctions(self):
    k_x, k_y, _ = self.cartesianKSymbols
    k_1, k_2, k_3 = self.grapheneKSymbols
    c_yz, c_zx, c_xy = self.orbitalShapeSymbols

    c_yz = ((1 + k_2 * k_3) / sp.sqrt(1 + k_x**2 + k_y**2)).subs(self.substitutionsKGraphene)
    c_zx = ((1 + k_1 * k_3) / sp.sqrt(1 + k_x**2 + k_y**2)).subs(self.substitutionsKGraphene)
    c_xy = ((1 + k_1 * k_2) / sp.sqrt(1 + k_x**2 + k_y**2)).subs(self.substitutionsKGraphene)
    return (c_yz, c_zx, c_xy)

  def __createCoordinateSystemAndDeltaTabs(self):
    self.N = CoordSys3D('N')
    self.deltaNearestTab = [
      0 * self.N.i - 1 * self.N.j,                                           # (0, -1)
      sp.sqrt(3)/sp.S(2) * self.N.i - sp.Rational(1, 2) * self.N.j,              # (√3/2, -1/2)
      sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(1, 2) * self.N.j,             # (√3/2, 1/2)
      -1 * (0 * self.N.i - 1 * self.N.j),                                    # (0, 1)
      -1 * (sp.sqrt(3)/sp.S(2) * self.N.i - sp.Rational(1, 2) * self.N.j),       # (-√3/2, 1/2)
      -1 * (sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(1, 2) * self.N.j)       # (-√3/2, -1/2)
    ]
    self.deltaNextTab = [
      sp.sqrt(3) * self.N.i + 0 * self.N.j,                                  # (√3, 0)
      sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(3, 2) * self.N.j,              # (√3/2, 3/2)
      -sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(3, 2) * self.N.j,             # (-√3/2, 3/2)
      -1 * (sp.sqrt(3) * self.N.i + 0 * self.N.j),                           # (-√3, 0)
      -1 * (sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(3, 2) * self.N.j),       # (-√3/2, -3/2)
      -1 * (-sp.sqrt(3)/sp.S(2) * self.N.i + sp.Rational(3, 2) * self.N.j)       # (√3/2, -3/2)
    ]

  def __getSingleBasisFunction(self, eigenvec):
    N = CoordSys3D('N')
    nOrbs = 3
    nSites = 6
    k_x, k_y, k_z = self.cartesianKSymbols
    kVec = self.N.i * k_x + self.N.j * k_y + self.N.k * k_z


    #Nearest neighbours
    orbBasisFunctions = []
    latexOrbBasisFunctions = []
    orbKxKyBasisFunctions = []
    for io in range(nOrbs):
      basisOrb = 0
      for ir in range(nSites):
        idx = io * nSites + ir
        basisOrb += sp.exp(-sp.I * (kVec & self.deltaNearestTab[ir])) * eigenvec[idx]
      trig = sp.simplify(basisOrb.rewrite(sp.cos))
      orbKxKyBasisFunctions.append(trig) # To solve for zeros later
      newBasis = trig.subs(self.substitutionsNearest)
      orbBasisFunctions.append(newBasis)
      latexOrbBasisFunctions.append(sp.latex(newBasis))
    display(orbBasisFunctions)
    print(latexOrbBasisFunctions)
    self.__getZerosOfBasisFunction(orbKxKyBasisFunctions)


    #Next nearest neighbors
    orbBasisFunctions = []
    latexOrbBasisFunctions = []
    orbKxKyBasisFunctions = []
    for io in range(nOrbs):
      basisOrb = 0
      for ir in range(nSites):
        idx = io * nSites + ir
        basisOrb += sp.exp(-sp.I * (kVec & self.deltaNextTab[ir])) * eigenvec[idx]
      trig = sp.simplify(basisOrb.rewrite(sp.cos))
      orbKxKyBasisFunctions.append(trig) # To solve for zeros later
      newBasis = trig.subs(self.substitutionsNext)
      orbBasisFunctions.append(newBasis)
      latexOrbBasisFunctions.append(sp.latex(newBasis))
    display(orbBasisFunctions)
    print(latexOrbBasisFunctions)
    self.__getZerosOfBasisFunction(orbKxKyBasisFunctions)


  def __getZerosOfBasisFunction(self, orbBasisFunctions):
    """
    Orbitals are stored in the following order:
    d_yz, d_zx, d_xy
    """
    k_x, k_y, _ = self.cartesianKSymbols
    c_yz, c_zx, c_xy = self.__getSymbolicOrbitalShapeFunctions()

    X, Y = self.kMesh

    equation = orbBasisFunctions[0]*c_yz + orbBasisFunctions[1]*c_zx + orbBasisFunctions[2]*c_xy
    display(equation)
    f = sp.lambdify((k_x, k_y), equation, "numpy")

    X, Y = self.kMesh

    Z = f(X, Y)
    Z = Z / np.max(np.abs(Z)) #Normalize

    for part in ("real", "imag"):
      zPlot = np.zeros(Z.shape)
      cbarLabel = ""
      if part == "real":
        zPlot = np.real(Z)
        cbarLabel = r"$\mathfrak{Re}\left( \Gamma  \right)$ (meV)"
      elif part == "imag":
        zPlot = np.imag(Z)
        cbarLabel = r"$\mathfrak{Im}\left( \Gamma  \right)$ (meV)"

      norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
      fig = plt.figure(figsize=(7, 5), dpi=100)
      # Set up GridSpec (1 row, 1 column, with some spacing)
      gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2, right=0.8, top=0.8, bottom=0.25)
      ax = fig.add_subplot(gs[0,0])
      cax = ax.pcolormesh(X, Y, zPlot, cmap='bwr', norm=norm)
      cbar = fig.colorbar(cax, ax=ax)
      cbar.set_label(cbarLabel)
      ax.contour(X, Y, zPlot, levels=[0], colors='black', linestyles='dotted', linewidths=0.5)

      self.__plotFirstBrillouinZoneBoundary()
      ax.set_xlabel(r"$k_x~(\tilde{a}^{-1})$")
      ax.set_ylabel(r"$k_y~(\tilde{a}^{-1})$")
      plt.show()
      plt.close()

  def __plotFirstBrillouinZoneBoundary(self):
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
    plt.plot(
      brillouinZoneVertices[:, 0],
      brillouinZoneVertices[:, 1],
      "--",
      color="black",
      linewidth=1,
    )

  def __plotRGBLegend(self):
    size = 512
    image = np.zeros((size, size, 3))

    # Triangle vertices for RGB channels
    v1 = np.array([1.0, 0.0])  # Red corner
    v2 = np.array([0.0, 0.0])  # Green corner
    v3 = np.array([0.5, np.sqrt(3)/2])  # Blue corner (top)

    # Coordinate grid
    x = np.linspace(0, 1, size)
    y = np.linspace(0, np.sqrt(3)/2, size)
    X, Y = np.meshgrid(x, y)

    # For each pixel, compute barycentric coordinates
    def barycentric_coords(x, y):
        # Transformation matrix for barycentric coordinates
        detT = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])
        l1 = ((v2[0] - x) * (v3[1] - y) - (v3[0] - x) * (v2[1] - y)) / detT
        l2 = ((v3[0] - x) * (v1[1] - y) - (v1[0] - x) * (v3[1] - y)) / detT
        l3 = 1.0 - l1 - l2
        return l1, l2, l3

    for i in range(size):
        for j in range(size):
            x_val, y_val = X[i, j], Y[i, j]
            l1, l2, l3 = barycentric_coords(x_val, y_val)
            if (l1 >= 0) and (l2 >= 0) and (l3 >= 0):  # Inside triangle
                image[i, j, :] = [l1, l2, l3]

    # Plot the triangle
    ax = plt.gca()
    ax.imshow(image, extent=(0, 1, 0, np.sqrt(3)/2), origin='lower')

    # Add triangle border and labels
    triangle = Polygon([v1, v2, v3], closed=True, edgecolor='k', fill=False, lw=2)
    ax.add_patch(triangle)
    ax.text(*v1, r'$d_{yz}$', color='black', ha='right', va='top')
    ax.text(*v2, r'$d_{zx}$', color='black', ha='left', va='top')
    ax.text(*v3, r'$d_{xy}$', color='black', ha='center', va='bottom')

    ax.axis('off')
    plt.tight_layout()
    plt.show()

  def __setPlotParams(self):
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