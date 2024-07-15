"Ce programme simule un problÃ¨me de convection 1D par Volumes Finis"
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
x = np.linspace(a,b,N)
u = np.zeros(N)

# Fonction pour la condition initiale
def g(x):
    x = x % b
    if (1/2)<= x <= (3/2):
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



# Nombre CFL (0 < CFL <= 1)
CFL = 0.8
# Calcul du pas de temps pour assurer la stabilité
dt = CFL * dx / abs(c)
lamda =dt / dx

unewp = np.zeros(N)
unewx = np.zeros(N)


plt.pause(0.5)
lamda = dt/dx
Tfinal = 2.5
temps = 0

while (temps < Tfinal):
    for i in range(N):
        unewx[i] = g(x[i] - c * temps)
    for i in range(1,N-1):
        flux_left = (c/2)*(u[i-1]+u[i])-1/2*abs(c)*(u[i]-u[i-1])
        flux_right=(c/2)*(u[i]+u[i+1])-1/2*abs(c)*(u[i+1]-u[i])
        unewp[i] = u[i] - lamda*(flux_right - flux_left)# SchÃ©ma dÃ©centrÃ© amont
    # Conditions aux limites de Neumann (DÃ©rivÃ©es nulles)
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
   
# Err = Uex - Un
Err = 0   
for i in range(0,N):
    Err = Err + dx*(ux[i]-u[i])**2
Err = np.sqrt(Err)
#print(Err)
# norm_l2 = norm(Err, 2)
#print(norm_l2)
print(N-1)
print(dx)
print(np.log(dx))
print(Err)
v = np.log(Err)
print(v)

