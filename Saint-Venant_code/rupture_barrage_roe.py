import numpy as np
import matplotlib.pyplot as plt

# Paramètres du domaine
a = 0
b = 1000
N = 300
dx = (b - a) / (N - 1)
x = np.linspace(a, b, N)

# Nombre CFL (0 < CFL <= 1)
CFL = 0.6

# Initialisation de la solution u et h
u = np.zeros(N)
h = np.zeros(N)
q=np.zeros(N)
qnewx=np.zeros(N)
qnewp=np.zeros(N)
# Fonction pour la condition initiale de la hauteur h (rupture du barrage)
def f(x):
    if 0 <= x <= 500:
        return 5  # Hauteur de l'eau avant la rupture
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
Tfinal = 30  # Temps final de la simulation
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
    P= np.ones((2, 2))
    D= np.zeros((2, 2))
    PI= np.ones((2, 2))
    y=np.zeros(2)
    
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
        h_tild_r=(h[i]+h[i+1])/2
        if h[i]==0 and h[i+1]==0:
            u_tild_r=0
        else:   
            u_tild_r=((np.sqrt(h[i])*u[i]+np.sqrt(h[i+1])*u[i+1])/(np.sqrt(h[i])+np.sqrt(h[i+1])))
        y[0]=h[i+1]-h[i]
        y[1]= h[i+1]*u[i+1]-h[i]*u[i]
        P[1][0]=u_tild_r-np.sqrt(g*h_tild_r)
        P[1][1]=u_tild_r+np.sqrt(g*h_tild_r)
        D[0][0]= np.abs(u_tild_r-np.sqrt(g*h_tild_r))
        D[1][1]= np.abs(u_tild_r+np.sqrt(g*h_tild_r))
        PI[0][0]=u_tild_r+np.sqrt(g*h_tild_r)
        PI[1][0]=np.sqrt(g*h_tild_r)-u_tild_r
        PI[0][1]=-1
        A = np.dot(P, np.dot(D, np.dot(PI, y)))
        if h_tild_r==0:
            flux_right_h=0
            flux_right_u =0
        else:
            flux_right_h = 1/2*(u[i]*h[i] +u[i+1]*h[i+1]) -(1/(4*np.sqrt(g*h_tild_r)))*A[0]
            flux_right_u = 1/2*((u[i]**2)*h[i] +(1/2)*g*h[i]**2 + (u[i+1]**2)*h[i+1] +(1/2)*g*h[i+1]**2)-(1/(4*np.sqrt(g*h_tild_r)))*A[1]
        h_tild_r=(h[i-1]+h[i])/2
        if h[i]==0 and h[i-1]==0:
            u_tild_r=0
        else:

            u_tild_r=((np.sqrt(h[i-1])*u[i-1]+np.sqrt(h[i])*u[i])/(np.sqrt(h[i-1])+np.sqrt(h[i])))
        y[0]=h[i]-h[i-1]
        y[1]= h[i]*u[i]-h[i-1]*u[i-1]
        P[1][0]=u_tild_r-np.sqrt(g*h_tild_r)
        P[1][1]=u_tild_r+np.sqrt(g*h_tild_r)
        D[0][0]= np.abs(u_tild_r-np.sqrt(g*h_tild_r))
        D[1][1]= np.abs(u_tild_r+np.sqrt(g*h_tild_r))
        PI[0][0]=u_tild_r+np.sqrt(g*h_tild_r)
        PI[1][0]=np.sqrt(g*h_tild_r)-u_tild_r
        PI[0][1]=-1
        A = np.dot(P, np.dot(D, np.dot(PI, y)))
        if h_tild_r==0:
            flux_left_h =0
            flux_left_u =0
        else:
            flux_left_h = 1/2*(u[i-1]*h[i-1] +u[i]*h[i]) -(1/(4*np.sqrt(g*h_tild_r)))*A[0]
            flux_left_u = 1/2*((u[i-1]**2)*h[i-1] +(1/2)*g*h[i-1]**2 + (u[i]**2)*h[i] +(1/2)*g*h[i]**2)-(1/(4*np.sqrt(g*h_tild_r)))*A[1]
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
          qnewx[i]=hnewx[i]*unewx[i]
          
      # Mise à jour du temps
    temps += dt
        
    u = unewp.copy()
    h = hnewp.copy()
    q=qnewp.copy()
    # Tracé des courbes de u et unewx au cours du temps
    plt.plot(x, q, '-r', label='Exact Solution' if temps <= dt else "")
    plt.plot(x,qnewx, '-b', label='Numerical Solution' if temps <= dt else "")
   # plt.plot(x, q, '-g', label='Numerical Solution' if temps <= dt else "")
    plt.grid()
    plt.pause(0.01)

plt.legend()
plt.show()

