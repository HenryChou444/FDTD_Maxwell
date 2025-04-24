import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Gaussian
A = 1.0  # Amplitude of the Gaussian
t_0 = 0  # Center of the Gaussian
sigma = 1.0  # Standard deviation (controls the width of the Gaussian)

# Create the time grid
t = np.linspace(-5, 5, 500)  # Time grid from -5 to 5 with 500 points

# Define the Gaussian function
gaussian = A * np.exp(-((t - t_0) ** 2) / (2 * sigma ** 2))

# Plot the Gaussian
plt.plot(t, gaussian, label=f"Gaussian (σ={sigma})")
plt.xlabel("t")
plt.ylabel("Amplitude")
plt.title("Gaussian Function")
plt.legend()
plt.grid()
plt.show()