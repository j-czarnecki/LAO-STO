import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern Roman'
plt.rcParams['font.sans-serif'] = 'Computer Modern Sans serif'
plt.rcParams['font.monospace'] = 'Computer Modern Typewriter'
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
# Optionally, add custom LaTeX preamble
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}'

b1 = [4 * np.pi/3* np.sqrt(3)/2, 4 * np.pi/3 * -0.5]
b2 = [0, 4 * np.pi/3]

# Plot the vectors with other shades of blue and red
plt.figure()
plt.quiver(0, 0, b1[0], b1[1], angles='xy', scale_units='xy', scale=1, color='red', label=r'$\vec{b}_1$')
plt.quiver(0, 0, b2[0], b2[1], angles='xy', scale_units='xy', scale=1, color='blue', label=r'$\vec{b}_2$')

#plt.quiver(b2[0], b2[1], b1[0] + b2[0], b1[1] + b2[0], angles='xy', scale_units='xy', scale=1, color='red', alpha = 1, linestyle = '--')
#plt.quiver(b1[0], b1[1], b2[0], b2[1], angles='xy', scale_units='xy', scale=1, color='blue', alpha = 1, linestyle = '--')

plt.text(1, 0.6, 'I BZ', fontsize = 30)

plt.fill_between([0, b1[0]], [0, b1[1]], [b1[0] + b2[0] + .48, b1[1] + b2[1]], color='purple', alpha=0.15)

# Set plot limits
plt.xlim(-5, 5)
plt.ylim(-5, 5)

# Add labels
plt.xlabel(r'$k_x (\tilde{a}^{-1})$', labelpad = 0)
plt.ylabel(r'$k_y (\tilde{a}^{-1})$')
#plt.grid()
plt.legend()

# Show plot
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.savefig('Basis.png')
#plt.show()