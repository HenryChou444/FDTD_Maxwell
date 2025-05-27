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
M = 300 # Number of x steps
N = 150 # Number of y steps
Q = 300  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
y = np.linspace(0, (M - 1) * dy, N)  # Not really used, but clearer that way
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
epsilon_r = np.ones((N, M))
#epsilon_r[150:199, 0:199] = e_r  # Set the permittivity of the dielectric region
sigma = np.zeros((N, M))  # Conductivity grid


# Interface parameters
angle_deg = 15  # Angle of the interface in degrees
angle_rad = np.deg2rad(angle_deg)
y_start = 20  # Start at the bottom


# Compute y_end so the interface spans the whole x axis
y_end = int(y_start + (M - 1) * np.tan(angle_rad))
y_end = min(y_end, N-1)  # Ensure we don't go out of bounds
print(f"y_end: {y_end}")

for n in range(N):
    if y_start <= n <= y_end:
        # The interface moves in x as y increases, forming a trapezoid
        x_interface = int((n - y_start) / np.tan(angle_rad))
        x_interface = min(x_interface, M-1)
        epsilon_r[n, x_interface:M] = e_r  # Set dielectric to the right of the interface
    elif n  < y_start:
        epsilon_r[n, :] = e_r  # All dielectric above the trapezoid


# Create Jz
omega = 2 * np.pi * f  # Angular frequency
J = np.zeros((N, M))  # Current density


# Create Fields
E = np.zeros((N, M))  # Electric field, last sample is (M-1)
E_max = np.zeros((N, M))  # Electric field amplitude
E_figure = np.zeros((N, M))  # Electric field figure
Bx = np.zeros((N-1, M))  #
By = np.zeros((N, M-1))  # Magnetic field, last sample is (M-2)


# Initialize the figure for 2D plotting
fig, ax = plt.subplots(figsize=(10, 8))
norm = Normalize(vmin=-1, vmax=1)  # Initial normalization (vmax will be updated dynamically)
im = ax.imshow(E, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm = norm)
cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.1, shrink=1, aspect=50)
cbar.set_label("$E_z$ [V/m]", fontsize=20)
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(0, (N - 1) * dy)
ax.set_xlabel("x [m]", fontsize=20)
ax.set_ylabel("y [m]", fontsize=20)
#ax.set_title("2D FDTD Simulation", fontsize=20)

# Create a mask for the dielectric zone
dielectric_mask = (epsilon_r == e_r).astype(float)
ax.imshow(
    dielectric_mask,
    extent=[0, (M - 1) * dx, 0, (N - 1) * dy],
    origin="lower",
    cmap="Greys",
    alpha=0.3,  # Adjust transparency as needed
    vmin=0, vmax=1
)

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
        for m in range(50, 200, 1):
            if frame < 39 :
                E[149, m] = np.sin(omega * frame * dt)  # Sine wave source
            else :
                E[149, m] = 0
        for n in range(1, N - 1):
            for m in range(1, M - 1): #1 compris, M-1 exclu
                if sigma[n, m] == -1: #use -1 instead of +inf
                    E[n, m] = 0
                elif sigma[n, m] == 0:
                    E[n, m] = E[n, m] + 1/a/epsilon_r[n,m] * (-Bx[n, m] + Bx[n-1, m] + By[n , m] - By[n, m-1]) - (J[n,m])/epsilon_r[n,m] #Normalized J
                    if abs(E[n, m]) > E_max[n, m]:
                        E_max[n, m] = abs(E[n, m])
                else:
                    E[n, m] = (1- sigma[n,m]/epsilon_r[n,m])*E[n, m] + 1/a/epsilon_r[n,m] * (-Bx[n, m] + Bx[n-1, m] + By[n , m] - By[n, m-1])
        print(f"E[N//2, M//2] at frame {frame}: {E[N//2, M//2]}")


        if frame < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
            for n in range(0, N - 1): 
                for m in range(1, M - 1): #0 compris, M-1 exclu
                    Bx[n, m] = Bx[n,m] - 1/a *(E[n+1, m] - E[n, m]) 

            for n in range(1, N - 1): 
                for m in range(0, M - 1): #0 compris, M-1 exclu
                    By[n, m] = By[n,m] + 1/a *(E[n, m + 1] - E[n, m])
        if frame == 225 : #save the image
            plt.savefig("2D_planar_wave_rectangle_domain.png", dpi=300, bbox_inches='tight')
    im.set_array(E)  # Update the color scale with the new E values
    return [im], ax.title

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=False, repeat=False)
# Save the animation
#ani.save("2D_planar_wave_rectangle_domain7.mp4", fps=45) #Must save before plt.show() but then additional waiting time
plt.show()
#print("Animation saved")
# Show the animation




# # Plot E_max in a 2D plot
# fig_max, ax_max = plt.subplots()

# # Normalize E_max to its maximum value
# norm_max = Normalize(vmin=np.percentile(E_max, 5), vmax= E_max[N//2, M//2 + 2])
# #norm_max = LogNorm(vmin=1e-5, vmax= np.max(E_max))
# # Create the 2D plot for E_max
# im_max = ax_max.imshow(E_max, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm = norm_max)
# #im_max = ax_max.imshow(E_max, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm=LogNorm(vmin=np.min(E_max[E_max > 0]), vmax=np.max(E_max)))
# cbar_max = fig_max.colorbar(im_max, ax=ax_max)
# cbar_max.set_label("Amplitude [V/m]", fontsize=14)

# # Set axis labels and title
# ax_max.set_xlim(0.14, 1.1)
# ax_max.set_ylim(0.14, 1.1)
# ax_max.set_xlabel("x [m]", fontsize=14)
# ax_max.set_ylabel("y [m]", fontsize=14)

# V = E_max * Lambda/np.pi # Voltage (V)
# P = (1/2) *(50/123**2) * V**2 # Power (W)
# # Normalize E_max to its maximum value
# norm_V = Normalize(vmin=0, vmax=V[N//2, M//2 + 2])  # Normalization for V
# norm_P = Normalize(vmin=0, vmax=P[N//2, M//2 + 2])  # Normalization for P

# # Plot V in a 2D plot
# fig_V, ax_V = plt.subplots()
# im_V = ax_V.imshow(V, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm = norm_V)
# #im_V = ax_V.imshow(V, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm=LogNorm(vmin=np.min(V[V > 0]), vmax=np.max(V)))
# cbar_V = fig_V.colorbar(im_V, ax=ax_V)
# cbar_V.set_label("Amplitude [V]", fontsize=20)

# # Set axis labels and title for V
# ax_V.set_xlim(0.14, 1.1)
# ax_V.set_ylim(0.14, 1.1)
# ax_V.set_xlabel("x [m]", fontsize=20)
# ax_V.set_ylabel("y [m]", fontsize=20)
# #ax_V.set_title("Voltage Distribution (V)", fontsize=14)

# # Plot P in a 2D plot
# fig_P, ax_P = plt.subplots()
# im_P = ax_P.imshow(P, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm=norm_P)
# #im_P = ax_P.imshow(P, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm=LogNorm(vmin=np.min(P[P > 0]), vmax=np.max(P)))
# cbar_P = fig_P.colorbar(im_P, ax=ax_P)
# cbar_P.set_label("Puissance [W]", fontsize=20)

# # Set axis labels and title for P
# ax_P.set_xlim(0.14, 1.1)
# ax_P.set_ylim(0.14, 1.1)
# ax_P.set_xlabel("x [m]", fontsize=20)
# ax_P.set_ylabel("y [m]", fontsize=20)
# #ax_P.set_title("Power Distribution (P)", fontsize=14)

# # Show both plots
# plt.show()
