import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import c
from scipy.constants import epsilon_0
from matplotlib.patches import Rectangle

# Paramètres
f = 2.4e9  # Fréquence (Hz)
Lambda = c / f  # Longueur d'onde (m)
dx = Lambda / 20  # Pas spatial (m)
a = 2 # constante qui permet de changer le pas temporel (non utilisée finalement)
dt = dx / (a * c)  # Pas temporel (s)
e_r = 4 # Permittivité relative du sol
M = 200 # Nombre de points spatiaux
Q = 800  # Nombre de pas de temps
x = np.linspace(0, (M - 1) * dx, M)  # Grille spatiale
t = np.linspace(0, (Q - 1) * dt, Q)  # Grille temporelle
# Grille de permittivité
wall_start = 3 * M // 5   # Début du mur
wall_width = M // 4   # Largeur du mur
epsilon_r = np.ones((1, M))
epsilon_r[0, wall_start : wall_start + wall_width] = e_r
sigma = np.zeros((1, M))  # Grille de conductivité
sigma[0, wall_start : wall_start + wall_width] = 0.1 # Conductivité du sol

# Création de Jz
omega = 2 * np.pi * f  # Pulsation angulaire
J = np.zeros((Q, M))  # Densité de courant
source_position = (M // 2)  # Position de la source (int)
J[:, source_position] = -np.sin(omega * t)  # Source sinusoïdale

# Création des champs
E = np.zeros((1, M))  # Champ électrique, dernier échantillon (M-1)
B = np.zeros((1, M-1))  # Champ magnétique, dernier échantillon (M-2)
E_boundary = np.zeros((2, 2))  # Conditions aux limites pour E

# Initialisation de la figure
fig, ax = plt.subplots()
animated_source, = ax.plot([], [])  # E
ax.set_xlim(0, (M - 1) * dx)
ax.set_ylim(-2, 2)
ax.set_xlabel("x [m]")
ax.set_ylabel("E [V/m]")
ax.set_title("FDTD 1D")

# Fonction de réinitialisation en fin d'animation
def reset():
    global E, B, E_boundary
    print("Réinitialisation de tous les vecteurs à zéro...")
    E[:] = 0  # Réinitialise le champ électrique
    B[:] = 0  # Réinitialise le champ magnétique
    E_boundary[:] = 0  # Réinitialise les conditions aux limites

# Fonction qui permet de calculer les champs et envoyer les données à l'animation
def update(frame):

    if frame == 0:  # Réinitialisation au début de l'animation
        reset()
    # Conditions aux limites
    if frame > 0:
        E_boundary[0, 0] = E_boundary[1, 0]  # Mise à jour de la condition limite gauche
        E_boundary[1, 0] = E[0, 1]  # Stocke la prochaine condition limite gauche
        E_boundary[0, 1] = E_boundary[1, 1]  # Mise à jour de la condition limite droite
        E_boundary[1, 1] = E[0, M-2]  # Stocke la prochaine condition limite droite
    # Calculs FDTD
    # NB : On ne sait pas faire de demi-indices pour B donc B[q + 1/2, m + 1/2] devient B[q, m]
    # Comme on ne stocke pas la matrice, q n'est pas utilisé, on met simplement à jour un vecteur ligne
        for m in range(1, M - 1): #1 compris, M-1 exclu
            if sigma[0, m] == 0:
                E[0, m] = E[0, m] + 1/a/epsilon_r[0,m] * (B[0, m]- B[0, m-1]) - (J[frame,m])/epsilon_r[0,m] #J est normalisé par epsilon_r
            else:
                E[0, m] = (1- sigma[0,m]/epsilon_r[0,m])*E[0, m] + 1/a/epsilon_r[0,m] * (B[0, m]- B[0, m-1])
   
    #   # Application des conditions aux limites
        if frame > a - 1 : # Passe la première itération
            E[0, 0] = E_boundary[0, 0] 
            E[0, M-1] =  E_boundary[0,1]

        if frame < Q-1: #Ne calcule pas la dernière ligne de B (pas utile pour trouver E)
            for m in range(0, M - 1): #0 compris, M-1 exclu
                B[0, m] = B[0,m] + 1/a *(E[0, m+1] - E[0, m]) 

   
    animated_source.set_data(x, E[0, :]) # Met à jour le champ E
    return animated_source, ax.title

# Create the animation
ani = FuncAnimation(fig, update, frames=Q, interval=15, blit=False, repeat=True)

# Save the animation

ani.save("1D_sine_source.mp4", fps=60)

# Show the animation
plt.show()
