import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
import os
# Get the directory of the current script
script_dir = os.path.dirname(__file__)

# Construct the full path to the file
file_path = os.path.join(script_dir, "E_matrix_1D_sine.npy")

# Load the matrix
E = np.load(file_path)


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

# Initialize the figure
fig, ax = plt.subplots()
animated_source, = ax.plot([], [])  # E
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("x [m]")
ax.set_ylabel("E [V/m]")
ax.set_title("1D FDTD")

# Animation function
#animated_source.set_data([x[source_position]], [J[frame, source_position]])  # Update source
def update(frame):
    animated_source.set_data(x, E[frame, :]) # Update E field
    return animated_source,

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=True)


# Show the animation
plt.show()