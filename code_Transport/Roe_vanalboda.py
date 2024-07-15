#Vanalboda

import numpy as np
import matplotlib.pyplot as plt

# Longueur du domaine
a = 0
b = 6
# Nombre de noeuds 
N = 101
# Le pas du maillage
dx = (b-a)/(N-1)
# Vitesse de transport
c = 2
x = np.linspace(a, b, N)
# On se donne un nombre CFL tq (0 < CFL <= 1)
CFL = 0.4
# Calcul du pas du temps de sorte à vérifier la condition de stabilité
dt = CFL * dx / abs(c)
u = np.zeros(N)

# Fonction limiteur de pente
def phm(x):
    return max(0,(x+x**2)/(1+x**2))

# Fonction pour la condition initiale
def g(x):
    x = x % b
    if 1/2 <= x <= 3/2:
        return 1
    else:
        return 0

# Initialisation de la condition initiale
for i in range(N):
    u[i] = g(x[i])

# Tracé de la condition initiale
plt.plot(x, u, '-b', label='Initial Condition')
plt.grid()
plt.legend()

# Calcul du pas de temps pour assurer la stabilité
dt = CFL * dx / abs(c)
unewp = np.zeros(N)
unewx = np.zeros(N)

plt.pause(0.5)
lamda = dt / dx
Tfinal = 3
temps = 0
r=np.zeros(N)
while temps < Tfinal:
    for i in range(N):
        unewx[i] = g(x[i] - c * temps)
    for i in range(1, N-1):
        if (u[i+1] - u[i]) != 0:
            r[i] = (u[i] - u[i-1]) / (u[i+1] - u[i])
            
        else:
            r[i] = 0
        if (u[i] - u[i-1]) != 0:
             r[i-1] = (u[i-1] - u[i-2]) / (u[i] - u[i-1])
             
        else:
             r[i-1] = 0
        if i==N-2:
            r[i+1]==0
        else:
           if (u[i+2] - u[i+1]) != 0:
               r[i+1] = (u[i+1] - u[i]) / (u[i+2] - u[i+1])
             
           else:
               r[i+1] = 0
        U_L_right = u[i] + 0.5 * phm(r[i]) * (u[i+1] - u[i])
        if i==N-2:
           U_R_right=u[N-1]
        else:
           U_R_right = u[i+1] - 0.5 * phm(r[i+1]) * (u[i+2] - u[i+1])
      
        U_L_left = u[i-1] + 0.5 * phm(r[i-1]) * (u[i] - u[i-1])
        U_R_left = u[i] - 0.5 * phm(r[i]) * (u[i+1] - u[i])
        flux_left = (c / 2) * (U_L_left + U_R_left) - 0.5 * abs(c) * (U_R_left - U_L_left)
        flux_right = (c / 2) * (U_L_right + U_R_right) - 0.5 * abs(c) * (U_R_right - U_L_right)
        unewp[i] = u[i] - lamda * (flux_right - flux_left)
    # Conditions aux limites de Neumann (Dérivées nulles)
    unewp[0] = u[N-1]
    unewp[N-1] = unewp[N-2]
    temps += dt
    u = unewp.copy()
    ux = unewx.copy()

    # Tracé des courbes de u et ux au cours du temps
    plt.plot(x, ux, '-r', label='Exact Solution' if temps <= dt else "")
    plt.plot(x, u, '-b', label='Numerical Solution' if temps <= dt else "")
    plt.grid()
    plt.pause(0.01)

plt.legend()
plt.show()

# Calcul de l'erreur
Err = 0   
for i in range(N):
    Err += dx * (ux[i] - u[i]) ** 2
Err = np.sqrt(Err)
print(N-1)
print(dx)
print(np.log(dx))
print(Err)
v = np.log(Err)
print(v)