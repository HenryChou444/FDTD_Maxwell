import numpy as np

e1 = 4
e2 = 1
theta_i = np.arcsin(np.sqrt(e2 / e1))
theta_i = np.degrees(theta_i)
print(f"theta_i: {theta_i}")