import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from tqdm import tqdm

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
Q = 400  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
y = np.linspace(0, (M - 1) * dy, N)  # Not really used, but clearer that way
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
epsilon_r = np.ones((N, M))
sigma = np.zeros((N, M))  # Conductivity grid


# Create Jz
J = np.zeros((N, M))  # Current density
omega = 2 * np.pi * f  # Angular frequency
theta = 30*180/np.pi #Angle of incidence
delta_phi = 2 * np.pi * (np.sin(theta) * dx / Lambda)
# Create Fields
E = np.zeros((N, M))  # Electric field, last sample is (M-1)
E_max = np.zeros((N, M))  # Electric field amplitude
Bx = np.zeros((N-1, M))  #
By = np.zeros((N, M-1))  # Magnetic field, last sample is (M-2)


# Initialize the figure for 2D plotting
fig, ax = plt.subplots()
norm = Normalize(vmin=-1, vmax=1)  # Initial normalization (vmax will be updated dynamically)
im = ax.imshow(E, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm = norm)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Intensité de $E_z$ [V/m]", fontsize=14)
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
        #Create planar wave source
        for m in range(10, 190, 1):
            if frame < 40 :
                E[199, m] = np.sin(omega * frame * dt)  # Sine wave source
            else :
                E[199, m] = 0
        # Compute fields
        for n in range(1, N - 1):
            for m in range(1, M - 1): #1 compris, M-1 exclu
                if sigma[0, m] == 0:
                    E[n, m] = E[n, m] + 1/a/epsilon_r[n,m] * (Bx[n, m]- Bx[n-1, m] + By[n , m] - By[n, m-1]) - (J[n,m])/epsilon_r[n,m] #Normalized J
                    # if abs(E[n, m]) > E_max[n, m]:
                    #     E_max[n, m] = abs(E[n, m])
                else:
                    E[0, m] = (1- sigma[0,m]/epsilon_r[0,m])*E[0, m] + 1/a/epsilon_r[0,m] * (B[0, m]- B[0, m-1])
                if abs(E[n, m]) > E_max[n, m]:
                    E_max[n, m] = abs(E[n, m])
        #print(f"E[M//2, M//2] at frame {frame}: {E[M//2, M//2]}")


        if frame < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
            for n in range(0, N - 1): 
                for m in range(1, M - 1): #0 compris, M-1 exclu
                    Bx[n, m] = Bx[n,m] + 1/a *(E[n+1, m] - E[n, m]) 

            for n in range(1, N - 1): 
                for m in range(0, M - 1): #0 compris, M-1 exclu
                    By[n, m] = By[n,m] + 1/a *(E[n, m + 1] - E[n, m])
    im.set_array(E)  # Update the color scale with the new E values
    return [im], ax.title

# Create the animation
ani = FuncAnimation(fig, update, frames = Q, interval=15, blit=False, repeat=False)
# Save the animation
ani.save("2D__planar_wave_double_front_zero_atfer.mp4", fps=45) #Must save before plt.show() but then additional waiting time
print("Animation saved")
# Show the animation
plt.show()
print(f"E_max: {np.max(E_max)}")


