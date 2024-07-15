import numpy as np
import matplotlib.pyplot as plt

# Longueur du domaine
a = -2
b = 2
# Nombre de noeuds
N = 100
# Le pas du maillage
dx = (b-a)/(N-1)

# Coordonnées x des noeuds
x = np.linspace(a, b, N)

# Nombre CFL (0 < CFL <= 1)
CFL = 0.3

# Initialisation de la solution u
u = np.zeros(N)
# Fonction limiteur de pente
def phm(x):
    return (x + abs(x))/(1+x)

# Fonction pour la condition initiale
def g(x):
    if x < -1:
        return 1
    elif -1 <= x <= 0:
        return -x
    else:
        return 0

# Initialisation de la condition initiale
for i in range(N):
    u[i] = g(x[i])

# Tracé de la condition initiale
plt.plot(x, u, '-b', label='Initial Condition')
plt.grid()
plt.legend()

# Tableau pour la nouvelle valeur de u
unewp = np.zeros(N)
# Tableau pour la solution exacte au temps courant
unewx = np.zeros(N)

# Temps final de la simulation
Tfinal = 4

# Temps initial
temps = 0
r = np.zeros(N)

while temps < Tfinal:
    # Calcul du pas de temps pour assurer la stabilité
    max_u = np.max(np.abs(u))
    if max_u == 0:
        max_u = 1e-6  # Éviter la division par zéro
    dt = CFL * dx / max_u
    lamda = dt / dx

    # Mise à jour de la solution exacte au temps courant
    if temps < 1:
        for i in range(N):
            if x[i] - temps < -1:
                unewx[i] = 1
            elif -1 <= x[i] <= 0:
                unewx[i] = -x[i] / (1 - temps)
            else:
                unewx[i] = 0
    else:
        for i in range(N):
            if (2 * x[i] + 1) < temps:
                unewx[i] = 1
            else:
                unewx[i] = 0

    # Mise à jour de la solution numérique
    for i in range(1, N-1):
        if (u[i+1] - u[i]) != 0:
            r[i] = (u[i] - u[i-1]) / (u[i+1] - u[i])
        else:
            r[i] = 0
        
        if (u[i] - u[i-1]) != 0:
            r[i-1] = (u[i-1] - u[i-2]) / (u[i] - u[i-1])
        else:
            r[i-1] = 0

        if i == N-2:
            r[i+1] = 0
        else:
            if (u[i+2] - u[i+1]) != 0:
                r[i+1] = (u[i+1] - u[i]) / (u[i+2] - u[i+1])
            else:
                r[i+1] = 0
        
        U_L_right = u[i] + 0.5 * phm(r[i]) * (u[i+1] - u[i])
        if i == N-2:
            U_R_right = u[N-1]
        else:
            U_R_right = u[i+1] - 0.5 * phm(r[i+1]) * (u[i+2] - u[i+1])
        
        U_L_left = u[i-1] + 0.5 * phm(r[i-1]) * (u[i] - u[i-1])
        U_R_left = u[i] - 0.5 * phm(r[i]) * (u[i+1] - u[i])
        
        S_RR = max(U_R_right, U_L_right)
        S_LR = min(U_R_right, U_L_right)
        
        if S_RR == S_LR:
                if S_LR > 0:
                    flux_right = (U_L_right**2) / 2
                elif  S_RR==0:
                    flux_right = 0
                else:
                    flux_right = (U_R_right **2) / 2
                
        else:
            if S_LR > 0:
                flux_right = (U_L_right ** 2) / 2
            elif S_LR <= 0 <= S_RR:
                flux_right = (S_RR * (U_L_right ** 2) / 2 - S_LR * (U_R_right ** 2) / 2 + S_LR * S_RR * (U_R_right - U_L_right)) / (S_RR - S_LR)
            else:
                flux_right = (U_R_right ** 2) / 2
        
        S_RL = max(U_R_left, U_L_left)
        S_LL = min(U_L_left, U_R_left)
        
        if S_RL == S_LL:
            if S_LL > 0:
                flux_left = (U_L_left**2) / 2
            elif  S_RL==0:
                flux_left = 0
            else:
                flux_left = (U_R_left**2) / 2
        else:
            if S_LL > 0:
                flux_left = (U_L_left ** 2) / 2
            elif S_LL <= 0 <= S_RL:
                flux_left = (S_RL * (U_L_left ** 2) / 2 - S_LL * (U_R_left ** 2) / 2 + S_LL * S_RL * (U_R_left - U_L_left)) / (S_RL - S_LL)
            else:
                flux_left = (U_R_left ** 2) / 2
        
        unewp[i] = u[i] - lamda * (flux_right - flux_left)

    # Conditions aux limites de Neumann (Dérivées nulles)
    unewp[N-1] = unewp[N-2]
    unewp[0] = unewp[1]

    # Mise à jour du temps
    temps += dt
    # Mise à jour de u
    u = unewp.copy()

    # Tracé des courbes de u et unewx au cours du temps
    plt.plot(x, unewx, '-r', label='Exact Solution' if temps <= dt else "")
    plt.plot(x, u, '-b', label='Numerical Solution' if temps <= dt else "")
    plt.grid()
    plt.pause(0.01)

plt.legend()
plt.show()

# Calcul de l'erreur
Err = 0
for i in range(N):
    Err += dx * (unewx[i] - u[i]) ** 2
Err = np.sqrt(Err)

print("Nombre de noeuds :", N-1)
print("Pas du maillage :", dx)
print("Log du pas du maillage :", np.log(dx))
print("Erreur :", Err)
print("Log de l'erreur :", np.log(Err))