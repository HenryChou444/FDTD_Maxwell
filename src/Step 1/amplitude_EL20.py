import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import os

# Get the directory of the current script
script_dir = os.path.dirname(__file__)

# Construct the full path for saving the file
file_path = os.path.join(script_dir, "E_matrix_L20.npy")

E = np.load(file_path)  # Load the electric field matrix
# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
M = 150 # Number of space steps
Q = 1500  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid

amplitude = np.zeros(M)  # Initialize amplitude array
for m in range(0, M):
    amplitude[m] = (np.max(E[int(Q/2) : , m]) - np.min(E[int(Q/2) :, m]))/2  
    #amplitude[m] = np.max(E[int(Q/2) : , m])

plt.figure()
plt.plot(x, amplitude)
plt.title("")
plt.xlabel("x [m]", fontsize=24)
plt.ylabel("Amplitude [V/m]", fontsize=24)
plt.xlim(0, (M - 1) * dx)
plt.ylim(0, 1.5)
# Increase the size of the tick labels
plt.tick_params(axis='both', which='major', labelsize=18)  # Adjust the font size as needed
plt.grid()
plt.show()