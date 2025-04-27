import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0

# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
M = 100 # Number of space steps
Q = 200  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
#print(f"dt/e0 = {dt/epsilon_0}")  # 1.1764677777163337

# Create Jz
omega = 2 * np.pi * f  # Angular frequency
J = np.zeros((Q, M))  # Current density
source_position = (M // 2)  # Position of the source (int)
J[:, source_position] = -np.sin(omega * t)  # Sine wave source

# Create Fields
E = np.zeros((1, M))  # Electric field, last sample is (M-1)
B = np.zeros((1, M-1))  # Magnetic field, last sample is (M-2)
E_boundary = np.zeros((2, 2))  # Boundary conditions for E field
# Initialize the figure
fig, ax = plt.subplots()
animated_source, = ax.plot([], [])  # E
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(-3, 3)
ax.set_xlabel("x [m]")
ax.set_ylabel("E [V/m]")
ax.set_title("1D FDTD")

# Reset function to set all vectors to zero
def reset():
    global E, B, E_boundary
    print("Resetting all vectors to zero...")
    E[:] = 0  # Reset the electric field
    B[:] = 0  # Reset the magnetic field
    E_boundary[:] = 0  # Reset the boundary conditions

# Animation function
def update(frame):

    if frame == 0:  # Reset at the start of the animation
        reset()
    
    if frame > 0:
        E_boundary[0, 0] = E_boundary[1, 0]  # Update left boundary
        E_boundary[1, 0] = E[0, 1]  # Store next left boundary
        E_boundary[0, 1] = E_boundary[1, 1]  # Update right boundary
        E_boundary[1, 1] = E[0, M-2]  # Store next right boundary
        for m in range(1, M - 1): #1 compris, M-1 exclu
            E[0, m] = E[0, m] + 1/a *(B[0, m]- B[0, m-1]) - (J[frame,m]) #Normalized J
        #print(f"E[1, :] at frame {frame}: {E[1, :]}")

    #   # Boundary conditions
        if frame > a - 1 : #Not mandatory, but avoids setting 0s
            E[0, 0] = E_boundary[0, 0] 
            E[0, M-1] =  E_boundary[0,1]

        if frame < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
            for m in range(0, M - 1): #0 compris, M-1 exclu
                B[0, m] = B[0,m] + 1/a *(E[0, m+1] - E[0, m]) 

    # Update the title with the current frame number   
    animated_source.set_data(x, E[0, :]) # Update E field
    ax.set_title(f"1D FDTD - Frame {frame}")
    #print(f"Updated frame {frame}")
    return animated_source, ax.title

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=False, repeat=True)

# Save the animation
#ani.save("1D_sine_source_lossless_dielectric.mp4", fps=60)

# Show the animation

# q_specific = 123  # Replace with the desired time step
# plt.figure()  # Create a new figure
# plt.plot(x, E[q_specific, :])
# plt.xlim(0, (M - 1) * dx)  # Same x-axis limits as the animation
# plt.ylim(-1.5, 1.5)  # Same y-axis limits as the animation
# plt.xlabel("x [m]")
# plt.ylabel("E [V/m]")
#plt.legend()
#plt.grid()


plt.show()
