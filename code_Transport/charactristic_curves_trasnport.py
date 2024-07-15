import matplotlib.pyplot as plt
import numpy as np

# Définir la fonction pour l'équation x = 2t + k
def equation(t, k):
    return 2 * t + k

# Définir l'intervalle de valeurs pour t
t_values = np.linspace(-10, 10, 11)

# Définir les différentes valeurs de k
k_values = np.linspace(-10, 20, 16)

# Tracer les courbes pour chaque valeur de k avec la couleur noire
for k in k_values:
    x_values = equation(t_values, k)
    plt.plot(x_values, t_values, color='black')

# Ajouter des légendes et des titres
plt.xlabel('x')
plt.ylabel('t')
plt.title('Characteristic Curves of the Transport Equation')
plt.grid(True)

# Enregistrer le graphique
plt.savefig('characteristic_curves_transport_equation.png')

# Afficher le graphique
plt.show()
