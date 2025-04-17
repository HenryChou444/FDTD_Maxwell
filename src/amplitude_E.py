import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c


E = np.load("E_matrix_1D_sine.npy")  # Load the electric field matrix
# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 50  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
M = 150 # Number of space steps
Q = 2000  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid

amplitude = np.zeros(M)  # Initialize amplitude array
for m in range(0, M):
    amplitude[m] = (np.max(E[int(Q/2) : , m]) - np.min(E[int(Q/2) :, m]))/2  
    #amplitude[m] = np.max(E[int(Q/2) : , m])

plt.figure()
plt.plot(x, amplitude)
plt.title("")
plt.xlabel("x [m]", fontsize=14)
plt.ylabel("Amplitude [V/m]", fontsize=14)
plt.xlim(0, (M - 1) * dx)
plt.ylim(0, 0.5)
plt.grid()
# Save the figure with high quality
plt.savefig("amplitude_plot.pdf", dpi=300, bbox_inches="tight")
plt.show()