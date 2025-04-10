import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0

# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
dt = dx / (2 * c)  # Time step (s)
M = 1000  # Number of space steps
Q = 1000  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
#print(f"dt/e0 = {dt/epsilon_0}")  # 1.1764677777163337

# Create Jz
omega = 2 * np.pi * f  # Angular frequency
J = np.zeros((Q, M))  # Current density
source_position = M // 2  # Position of the source (int)
J[:, source_position] = np.sin(omega * t)  # Sine wave source

# Initialize the figure
fig, ax = plt.subplots()
animated_source, = ax.plot([], [], 'o', markersize=4, color='red')  # Sine wave source
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("x [m]")
ax.set_ylabel("E [V/m]")
ax.set_title("1D FDTD")

# Animation function
def update(frame):
    animated_source.set_data([x[source_position]], [J[frame, source_position]])  # Update source
    return animated_source,

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=5, blit=True)

# Show the animation
plt.show()