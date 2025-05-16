import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# Get screen size in pixels (works on many systems)
try:
    import tkinter as tk
    root = tk.Tk()
    width_px = root.winfo_screenwidth()
    height_px = root.winfo_screenheight()
    root.destroy()
except:
    width_px, height_px = 1920, 1080  # fallback
dpi = 100
# Parameters
f = 2.4e9  # Frequency (Hz)
Lambda = c / f  # Wavelength (m)
dx = Lambda / 20  # Spatial step (m)
dy = Lambda / 20  # Spatial step (m)
a = 2
dt = dx / (a * c)  # Time step (s)
e_r = 4 # Relative permittivity of ground
M = 400 # Number of x steps
N = 20 # Number of y steps
Q = 500  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
y = np.linspace(0, (M - 1) * dy, N)  # Not really used, but clearer that way
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
epsilon_r = np.ones((N, M))
sigma = np.zeros((N, M))  # Conductivity grid

# Add perfect conductor wall
sigma[N//2 + 6, :] = -1  # -1 instead of +inf for the code
sigma[N//2 - 6, :] = -1  # Perfect conductor wall

# Create Jz
omega = 2 * np.pi * f  # Angular frequency
J = np.zeros((N, M))  # Current density
source_position_x = (M // 2)  # Position of the source (int)
source_position_y = (N // 2)  # Position of the source (int)


# Create Fields
E = np.zeros((N, M))  # Electric field, last sample is (M-1)
E_max = np.zeros((N, M))  # Electric field amplitude
Bx = np.zeros((N-1, M))  #
By = np.zeros((N, M-1))  # Magnetic field, last sample is (M-2)


# Initialize the figure for 2D plotting
fig, ax = plt.subplots(figsize=(width_px/dpi, height_px/dpi), dpi=dpi)
norm_max = Normalize(vmin=-0.4682158550076026, vmax= 0.4682158550076026)
im = ax.imshow(E, extent=[0, (M - 1) * dx, 0, (N - 1) * dy], origin="lower", cmap="jet", norm = norm_max)
cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, shrink=1, aspect=100)
cbar.set_label("$E_z$ [V/m]", fontsize=25)
cbar.ax.tick_params(labelsize=20)
ax.set_xlim(0, (M - 1) * dx)
ax.set_xlabel("x [m]", fontsize=25)
ax.tick_params(axis='x', labelsize=20)
#ax.set_ylabel("y [m]", fontsize=25)
ax.set_yticklabels([])
ax.set_title("2D FDTD Simulation", fontsize=25)

# Center the plot and remove extra white space
fig.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15)

# Draw horizontal red lines for the conductor walls 
y1 = (N // 2 + 6) * dy
y2 = (N // 2 - 6) * dy
ax.set_ylim(y2, y1)
# line1 = ax.axhline(y=y1, color='red', linewidth=4, label='Mur Conducteur')
# line2 = ax.axhline(y=y2, color='red', linewidth=5)
#ax.legend(loc='upper right')

#1D Plot
fig_cut_anim, ax_cut_anim = plt.subplots()
line_cut, = ax_cut_anim.plot([], [], color="blue")
ax_cut_anim.set_xlim(0, (M - 1) * dx)
ax_cut_anim.set_ylim(-0.7, 0.7)  # You may want to update this after the first run
ax_cut_anim.set_xlabel("x [m]", fontsize=20)
ax_cut_anim.set_ylabel("$E_z$ [V/m]", fontsize=20)
ax_cut_anim.set_title("Coupe 1D entre les 2 plaques", fontsize=20)
ax_cut_anim.grid()

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
        J[source_position_y, source_position_x] = -np.sin(omega * frame * dt) # Sine wave source
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
    # if frame == 408 : #save the image 
    #     plt.savefig("2D_TE_1_norm_test.png", dpi=300, bbox_inches='tight')

    im.set_array(E)  # Update the color scale with the new E values
    E_cut = E[N // 2, :]
    line_cut.set_data(x, E_cut)
    return [im, line_cut], ax.title


# Create the animation
#ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=False, repeat=False)
ani = FuncAnimation(fig_cut_anim, update, frames=Q, interval=15, blit=False, repeat=False)
# Save the animation
#ani.save("2D_TE_1_norm_test.mp4", fps=45) #Must save before plt.show() but then additional waiting time
# Print min and max of E_max
plt.show()
# print("E_max min:", np.min(E_max))
# print("E_max max:", np.max(E_max))

