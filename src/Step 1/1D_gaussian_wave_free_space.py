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
Q = 800  # Number of time steps
x = np.linspace(0, (M - 1) * dx, M)  # Space grid
t = np.linspace(0, (Q - 1) * dt, Q)  # Time grid
#print(f"dt/e0 = {dt/epsilon_0}")  # 1.1764677777163337


# Create Jz
J = np.zeros((Q, M))  # Current density
source_position = M // 2  # Position of the source (int)

# Parameters for the Gaussian waveform
sigma = 10*dt # Standard deviation (controls the width of the Gaussian)
t_0 = 3*sigma # Center of the Gaussian 
J[:, source_position] = -np.exp(-((t - t_0) ** 2) / (2 * sigma ** 2)) # Gaussian source


# Create E field
E = np.zeros((Q, M))  # Electric field, last sample is (Q-1, M-1)
B = np.zeros((Q-1, M-1))  # Magnetic field, last sample is (Q-2, M-2)

for q in range(1, Q): # B_(q'+1/2) [m' + 1/2] = B_q [m] ; début à q=1, car C.I nulles

#   # Boundary conditions
    # if q > a - 1 :
    #     E[q, 0] = E[q-a, 1] 
    #     E[q, M-1] =  E[q-a, M-2]
   
    for m in range(1, M - 1): #1 compris, M-1 exclu
        #E[q, m] = E[q - 1, m] + 1/a *(B[q-1, m]- B[q-1, m-1]) - (dt / epsilon_0) * (J[q-1,m])
        E[q, m] = E[q - 1, m] + 1/a *(B[q-1, m]- B[q-1, m-1]) - (J[q-1,m]) #Normalized J

    if q < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
        for m in range(0, M - 1): #0 compris, M-1 exclu
            B[q, m] = B[q - 1,m] + 1/a *(E[q, m+1] - E[q, m])   



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

# Save the animation
#ani.save("1D_gaussian_source_free_spac_boundary_zero.mp4", fps=60)

q_specific = 80  # Replace with the desired time step
plt.figure()  # Create a new figure
plt.plot(x, E[q_specific, :], linewidth=3)
plt.xlim(0, (M - 1) * dx)  # Same x-axis limits as the animation
plt.ylim(-0.5, 1.5)  # Same y-axis limits as the animation
plt.xlabel("x [m]",fontsize=30)
plt.ylabel("E [V/m]", fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=25)  # Adjust the font size as needed

# Show the animation
plt.show()
