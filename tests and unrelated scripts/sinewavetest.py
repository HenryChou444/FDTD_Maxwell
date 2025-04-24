import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c

# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
dt = dx / (2 * c)  # Time step (s)
M = 1000  # Number of space steps
Q = 200  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid

# Create the sine wave
k = 2 * np.pi / Lambda  # Wave number
omega = 2 * np.pi * f  # Angular frequency

# Initialize the figure
fig, ax = plt.subplots()
line, = ax.plot(x, np.zeros_like(x))  # Line to animate
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("Space (m)")
ax.set_ylabel("Amplitude")
ax.set_title("1D Sine Wave Animation")

# Animation function
def update(frame):
    y = np.sin(k * x - omega * t[frame])  # Sine wave at time t[frame]
    line.set_ydata(y)
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=50, blit=True)

# Show the animation
plt.show()