import numpy as np 
import matplotlib.pyplot as plt


H_real = np.loadtxt('../OutputData/H_real.dat').view(float)
H_imag = np.loadtxt('../OutputData/H_imag.dat').view(float)
plt.matshow(H_real, cmap = 'seismic')
plt.colorbar()
plt.grid(1)
plt.savefig("../Plots/H_real.png")
plt.matshow(H_imag, cmap = 'seismic')
plt.colorbar()
plt.grid(1)
plt.savefig("../Plots/H_imag.png")
