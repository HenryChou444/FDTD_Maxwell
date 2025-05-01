import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from matplotlib.patches import Rectangle

# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
dy = Lambda / 20  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
e_r = 4 # Relative permittivity of ground
M = 200 # Number of x steps
N = 200 # Number of y steps
Q = 300  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
y = np.linspace(0, (M - 1) * dy, N)  # Not really used, but clearer that way
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
# Permittivity grid
#wall_start = 3 * M // 5   # Start of the wall
#wall_width = M // 4   # Width of the wall
epsilon_r = np.ones((N, M))
#epsilon_r[0, wall_start : wall_start + wall_width] = e_r
sigma = np.zeros((N, M))  # Conductivity grid
#sigma[0, wall_start : wall_start + wall_width] = 0.1 # Conductivity of the ground
#print(f"dt/e0 = {dt/epsilon_0}")  # 1.1764677777163337

# Create Jz
omega = 2 * np.pi * f  # Angular frequency
J = np.zeros((N, M))  # Current density
source_position = (M // 2)  # Position of the source (int)


# Create Fields
E = np.zeros((N, M))  # Electric field, last sample is (M-1)
Bx = np.zeros((N-1, M))  #
By = np.zeros((N, M-1))  # Magnetic field, last sample is (M-2)


# Initialize the figure for 2D plotting
fig, ax = plt.subplots()
im = ax.imshow(E, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="viridis")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Amplitude de $E_z$ [V/m]", fontsize=14)
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(0, (N - 1) * dy)
ax.set_xlabel("x [m]", fontsize=14)
ax.set_ylabel("y [m]", fontsize=14)
ax.set_title("2D FDTD Simulation", fontsize=14)
# ax.set_title("2D FDTD")
# ax.add_patch(wall)
#wall_height = 4  # Height of the wall (adjust as needed)
# wall = Rectangle(((wall_start-1) * dx , -wall_height), wall_width*dx, 2 * wall_height, linewidth=1, edgecolor='red', facecolor='none', linestyle='--', label="Wall")
# ax.add_patch(wall)
#ax.legend()

# Reset function to set all vectors to zero
def reset():
    global E, B, E_boundary
    print("Resetting all vectors to zero...")
    E[:] = 0  # Reset the electric field
    Bx[:] = 0  # Reset the magnetic field
    By[:] = 0  # Reset the magnetic field


# Animation function
def update(frame):

    if frame == 0:  # Reset at the start of the animation
        reset()
    if frame > 0:
        J[source_position, source_position] = -np.sin(omega * frame * dt)  # Sine wave source
        for n in range(1, N - 1):
            for m in range(1, M - 1): #1 compris, M-1 exclu
                if sigma[0, m] == 0:
                    E[n, m] = E[n, m] + 1/a/epsilon_r[n,m] * (Bx[n, m]- Bx[n-1, m] + By[n , m] - By[n, m-1]) - (J[n,m])/epsilon_r[n,m] #Normalized J
                else:
                    E[0, m] = (1- sigma[0,m]/epsilon_r[0,m])*E[0, m] + 1/a/epsilon_r[0,m] * (B[0, m]- B[0, m-1])
        print(f"E[M//2, M//2] at frame {frame}: {E[M//2, M//2]}")


        if frame < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
            for n in range(0, N - 1): 
                for m in range(1, M - 1): #0 compris, M-1 exclu
                    Bx[n, m] = Bx[n,m] + 1/a *(E[n+1, m] - E[n, m]) 

            for n in range(1, N - 1): 
                for m in range(0, M - 1): #0 compris, M-1 exclu
                    By[n, m] = By[n,m] + 1/a *(E[n, m + 1] - E[n, m])

    # Update the title with the current frame number   
    
    #ax.set_title(f"1D FDTD - Frame {frame}")
    #print(f"Updated frame {frame}")

    # # Figure at a specific time step
    # if frame == 605 :  # Replace with the desired time step
    #     # Create a new figure
    #     fig_new, ax_new = plt.subplots()
    #     ax_new.plot(x, E[0, :], label=f"Frame {frame}")
    #     ax_new.set_xlim(0, (M - 1) * dx)  # Same x-axis limits as the animation
    #     ax_new.set_ylim(-1.5, 1.5)  # Same y-axis limits as the animation
    #     ax_new.set_xlabel("x [m]", fontsize=24)
    #     ax_new.set_ylabel("E [V/m]", fontsize=24)
    #     ax_new.tick_params(axis='both', which='major', labelsize=14)

    #     # Add the rectangle patch for the wall
    #     wall_new = Rectangle(((wall_start - 1) * dx, -wall_height), wall_width * dx, 2 * wall_height,
    #                         linewidth=1, edgecolor='red', facecolor='none', linestyle='--', label="Wall")
    #     ax_new.add_patch(wall_new)

    #     # Show the new plot
    #     plt.show()µ
     # Update the plot with the magnitude of E
    im.set_array(E)  # Update the color scale with the new E values
    return [im], ax.title

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=False, repeat=True)

# Save the animation
ani.save("2D_sine_source_free_space.mp4", fps=45) #Must save before plt.show() but then additional waiting time
print("Animation saved")
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


