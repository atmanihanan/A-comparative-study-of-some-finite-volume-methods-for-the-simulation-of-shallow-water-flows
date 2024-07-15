import numpy as np
import matplotlib.pyplot as plt

# Paramètres du domaine
a = 0
b = 1000
N = 300
dx = (b - a) / (N - 1)
x = np.linspace(a, b, N)

# Nombre CFL (0 < CFL <= 1)
CFL = 0.5

# Initialisation de la solution u et h
u = np.zeros(N)
h = np.zeros(N)
q=np.zeros(N)
qnewp=np.zeros(N)
# Fonction pour la condition initiale de la hauteur h (rupture du barrage)
def f(x):
    if 0 <= x <= 500:
        return 5 # Hauteur de l'eau avant la rupture
    else:
        return 0  # Hauteur de l'eau après la rupture

# Initialisation de la condition initiale
for i in range(N):
    h[i] = f(x[i])

# Tracé de la condition initiale
plt.plot(x, h, '-r', label='Hauteur initiale (h)')
plt.xlabel('x')
plt.ylabel('Hauteur de l\'eau')
plt.title('Condition initiale de la rupture du barrage')
plt.grid()
plt.legend()
plt.show()

# Tableau pour les nouvelles valeurs de u et h
unewp = np.zeros(N)
hnewp = np.zeros(N)
unewx=np.zeros(N)
hnewx=np.zeros(N)
g = 9.81  # Accélération due à la gravité
Tfinal = 20  # Temps final de la simulation
temps = 0
while temps < Tfinal:
    # Calcul du pas de temps pour assurer la stabilité
    s_1 = np.abs(u + np.sqrt(g * h))
    s_2 = np.abs(u - np.sqrt(g * h))
    S_1max=max(s_1)
    S_2max=max(s_2)
    max_u = max(S_1max, S_2max)
   
    dt = CFL * dx / max_u
    lamda = dt / dx
    #mise a jour de solution exacte 
    for i in range(N):
        if x[i]<=(500-temps*np.sqrt(g*5)):
            hnewx[i]= 5
            unewx[i]= 0
            
        elif (500-temps*np.sqrt(g*5))<x[i]<(500+2*temps*np.sqrt(g*5)):
            unewx[i]= (2/3)*(np.sqrt(g*5)+(x[i]-500)/temps)
            hnewx[i]= (1/(9*g))*((2*np.sqrt(g*5)-(x[i]-500)/temps)**2)
        else:
            hnewx[i]=0
            unewx[i]=0
    # Mise à jour de la solution numérique
    for i in range(1, N-1):

        S_RR = max(u[i+1] +np.sqrt(g*h[i+1]),u[i] +np.sqrt(g*h[i]))
        S_LR = min(u[i+1] -np.sqrt(g*h[i+1]),u[i] -np.sqrt(g*h[i]))
        if S_RR == S_LR:
            if S_LR > 0:
                flux_right_h = u[i]*h[i]
                flux_right_u = (u[i]**2)*h[i] +(1/2)*g*h[i]**2
            elif  S_RR==0:
                flux_right_h = 0
                flux_right_u = 0
                
            else:
                flux_right_h = u[i+1]*h[i+1]
                flux_right_u = (u[i+1]**2)*h[i+1] +(1/2)*g*h[i+1]**2
            
        else:
            if S_LR > 0:
                flux_right_h = u[i]*h[i]
                flux_right_u = (u[i]**2)*h[i] +(1/2)*g*h[i]**2
            elif S_LR <= 0 <= S_RR:
                flux_right_h = (S_RR * u[i]*h[i] - S_LR * u[i+1]*h[i+1] + S_LR * S_RR * ( h[i+1]- h[i] )) / (S_RR - S_LR)
                flux_right_u =( S_RR * ((u[i]**2)*h[i] +(1/2)*g*h[i]**2) - S_LR * ((u[i+1]**2)*h[i+1] +(1/2)*g*h[i+1]**2) + S_LR * S_RR * ( u[i+1]- u[i] )) / (S_RR - S_LR)
            else:
                flux_right_h = u[i+1]*h[i+1]
                flux_right_u =(u[i+1]**2)*h[i+1] +(1/2)*g*h[i+1]**2
                
        S_RL = max(u[i] +np.sqrt(g*h[i]),u[i-1] +np.sqrt(g*h[i-1]))
        S_LL = min(u[i-1] -np.sqrt(h[i-1]),u[i] -np.sqrt(g*h[i]))
        if S_RL == S_LL:
            if S_LL > 0:
                flux_left_h = u[i-1]*h[i-1]
                flux_left_u =(u[i-1]**2)*h[i-1] +(1/2)*g*h[i-1]**2
            elif  S_RL==0:
                flux_left_h = 0
                flux_left_u= 0
            else:
                flux_left_h= u[i]*h[i]
                flux_left_u= (u[i]**2)*h[i] +(1/2)*g*h[i]**2
        else:
            if S_LL > 0:
                flux_left_h =u[i-1]*h[i-1]
                flux_left_u = (u[i-1]**2)*h[i-1] +(1/2)*g*h[i-1]**2
            elif S_LL <= 0 <= S_RL:
                flux_left_h =(S_RL * u[i-1]*h[i-1] - S_LL * u[i]*h[i] + S_LL * S_RL * ( h[i]- h[i-1] )) / (S_RL - S_LL)
                flux_left_u =  (S_RL * ((u[i-1]**2)*h[i-1] +(1/2)*g*h[i-1]**2) - S_LL * ((u[i]**2)*h[i] +(1/2)*g*h[i]**2)+ S_LL * S_RL * ( u[i]- u[i-1] )) / (S_RL - S_LL)
            else:
                flux_left_h = u[i]*h[i]
                flux_left_u = (u[i]**2)*h[i] +(1/2)*g*h[i]**2
        hnewp[i] = h[i] - lamda * (flux_right_h - flux_left_h)
        if hnewp[i]==0:
            unewp[i]=0
        else:
            unewp[i] = (1/hnewp[i])*(h[i]*u[i] - lamda * (flux_right_u - flux_left_u))

    # Conditions aux limites de Neumann (Dérivées nulles)
    hnewp[N-1] = hnewp[N-2]
    hnewp[0] = hnewp[1]
    unewp[N-1] = unewp[N-2]
    unewp[0] = unewp[1]
    for i in range(N):
          qnewp[i]=hnewp[i]*unewp[i]
          
      # Mise à jour du temps
    temps += dt
      # Mise à jour de u et h
    u = unewp.copy()
    h = hnewp.copy()
    q=qnewp.copy()
    # Tracé des courbes de u et unewx au cours du temps
    plt.plot(x, h, '-r', label='Exact Solution' if temps <= dt else "")
    #plt.plot(x, u, '-b', label='Numerical Solution' if temps <= dt else "")
    plt.plot(x, hnewx, '-g', label='Numerical Solution' if temps <= dt else "")
    plt.grid()
    plt.pause(0.01)

plt.legend()
plt.show()


