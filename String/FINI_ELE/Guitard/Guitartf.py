import numpy as np
import matplotlib.pyplot as plt

# Paramètres de la corde
L = 0.5               # Longueur de la corde (m)
c = 2600

rho_volume = 1150   #kg·m-3           
diameter = 0.00055 #m 
# The string is cicular
rho = rho_volume*(np.pi * diameter * diameter )/ 4   # Mass linéique kg.m
T = c**2*rho          # Tension de la corde (N)
print("the tension",T)
n_elements = 100       # Nombre d'éléments finis
dt = 10**(-8)             # Pas de temps (s)
t_end = 0.01          # Temps de simulation (s)
n_steps = int(t_end / dt)

# Discrétisation
n_nodes = n_elements + 1
dx = L / n_elements    # Longueur d'un élément
x = np.linspace(0, L, n_nodes)

# Matrices locales pour un élément
K_local = (T / dx) * np.array([[1, -1], [-1, 1]])  # Matrice de raideur locale
M_local = (rho * dx / 6) * np.array([[2, 1], [1, 2]])  # Matrice de masse locale

# Assemblage global des matrices K et M
K = np.zeros((n_nodes, n_nodes))

M = np.zeros((n_nodes, n_nodes))
for i in range(n_elements):
    dofs = [i, i + 1]
    K[np.ix_(dofs, dofs)] += K_local
    M[np.ix_(dofs, dofs)] += M_local
    
# Conditions aux limites : nœuds fixes
# ------------Modif K : représente les relations entre les déplacements nodaux et les forces internes dans la corde
K[0, :] = 0 #met à zéro la première ligne, annulant toute influence du déplacement du premier nœud u0
K[-1, :] = 0 #met à zéro la dernière ligne, annulant toute influence du déplacement du dernier nœud
# Fixer la diagonale à 1 pour ces nœuds 
K[0, 0] = 1 #Cela rend les déplacements de ces nœuds fixes égaux à zéro en garantissant une contrainte forte dans le système global
K[-1, -1] = 1
# ------------Modification de la matrice de masse globale MM
M[0, :] = 0 # Cela élimine la contribution des nœuds fixes à la dynamique du système.
M[-1, :] = 0
# Cela empêche toute contribution des forces d'inertie à ces nœuds fixes.
M[0, 0] = 1
M[-1, -1] = 1
step_mov = 0.003
#%% To change
# Initialisation des vecteurs
def diffQui(n_nodes,n_steps, step_mov,K,M, dt):  
    simulation = np.zeros((n_steps, n_nodes))      # Déplacements (position)
    # Ajouter une force initiale (par exemple, un impact au centre de la corde)

    # Boucle temporelle (méthode explicite de Newmark ou central differences)
    
    F = np.zeros(n_nodes)
    
    for time in range(0,n_steps-1):
        if time == 0 :
            F[n_nodes//4] = step_mov
        # Calcul des accélérations (F = M * a = K * U + F_ext)
        A = np.linalg.inv(M) @ (F - K @ simulation[time])
        
        # Mise à jour des positions (méthode central difference)
        simulation[time + 1] = 2 * simulation[time] - simulation[time - 1] + (dt**2) * A
        F = np.zeros(n_nodes)
    return simulation 

# Visualisation de la propagation des ondes
sim = diffQui(n_nodes,n_steps,step_mov,K,M, dt)
X = np.linspace(0, L, int(L/n_nodes))
 

time_steps_to_plot = range(0, n_steps, max(1, n_steps // 30))  # Select up to 30 frames to plot
for t in time_steps_to_plot:
    plt.plot(X, sim[:, t])
plt.ylabel('hauteur en m') # Légendes en x et y
plt.xlabel('x en m')
plt.title('hauteur de la corde en fonction de x pour différents temps')
plt.show()
