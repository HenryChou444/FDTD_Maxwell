import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
dy = Lambda / 20  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
e_r = 4 # Relative permittivity of ground
M = 400 # Number of x steps
N = 22 # Number of y steps
Q = 500  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
y = np.linspace(0, (M - 1) * dy, N)  # Not really used, but clearer that way
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
epsilon_r = np.ones((N, M))
sigma = np.zeros((N, M))  # Conductivity grid
b_2 = 10
# Add perfect conductor wall
sigma[N//2 + b_2, :] = -1  # -1 instead of +inf for the code
sigma[N//2 - b_2, :] = -1  # Perfect conductor wall

# things to compute 
beta = 2 * np.pi / Lambda
fc2 = c / (2*2*b_2*dy)
Lambda_g2 = 2*np.pi / (beta*np.sqrt(1 - (fc2/f)**2))

print(f"beta = {beta:.4e} [m]")
print(f"fc2 = {fc2:.4e} [Hz]")
print(f"Lambda_g2 = {Lambda_g2:.4e} [m]")