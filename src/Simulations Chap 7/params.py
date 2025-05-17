import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from scipy.constants import mu_0
from matplotlib.patches import Rectangle

# Parameters
f = 2.4e9  # Frequency (Hz)
w = 2 * np.pi * f  # Angular frequency (rad/s)
Lambda = c / f  # Wavelength (m)
e_r = 4 # Relative permittivity of ground
sigma = 0.1 # Conductivity of the ground
beta = w * np.sqrt(e_r * epsilon_0 * mu_0 / 2) * np.sqrt(np.sqrt(1 + (sigma / (w * e_r * epsilon_0))**2)  + 1 )  # Propagation constant
Lambda_eff = 2 * np.pi / beta  # Effective wavelength (m)
print(f"Lambda_eff = {Lambda_eff}")  # Effective wavelength