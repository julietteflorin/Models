import numpy as np
import matplotlib.pyplot as plt
import time
import wave
import math 


# Perturbation initiale (forme triangulaire)
# Déplacement initial qui ressemble à une corde de guitare pincée
# Perturbation initiale (forme triangulaire) avec contrôle de la hauteur
def perturbation_initiale(x, L, position_pincement, hauteur_max):
    # position_pincement : Position relative où la corde est pincée (par exemple, 0.25 pour un quart de la longueur)
    # hauteur_max : Hauteur maximale du pincement
    perturbation = np.zeros_like(x)  # Convertir en index
    for i in range(len(x)):
        if i <= position_pincement:  # Pente ascendante
            perturbation[i] = hauteur_max * (i / position_pincement)
        else:  # Pente descendante
            perturbation[i] = hauteur_max * ((len(x) - i) / (len(x) - position_pincement))
    return perturbation


def FEMgui(L, rho, T,c, n_elements, dx, n_nodes, position_pincement, hauteur_max, t_end, dt, n_steps):
    pression_in_time = np.zeros(n_steps)
    x = np.linspace(0, L, n_nodes)
    print(n_steps)
    dt_max = dx / c        # Condition CFL
    dt = min(dt, dt_max)
    print(f"dt utilisé : {dt}, dt_max (CFL) : {dt_max}")

    # Matrices locales
    K_local = (T / dx) * np.array([[1, -1], [-1, 1]])  # Matrice de raideur locale
    M_local = (rho * dx / 6) * np.array([[2, 1], [1, 2]])  # Matrice de masse locale

    # Assemblage global des matrices K et M
    K = np.zeros((n_nodes, n_nodes))
    M = np.zeros((n_nodes, n_nodes))
    for i in range(n_elements):
        dofs = [i, i + 1]
        K[np.ix_(dofs, dofs)] += K_local
        M[np.ix_(dofs, dofs)] += M_local

    # Conditions aux limites
    K[0, :] = 0
    K[-1, :] = 0
    K[0, 0] = 1
    K[-1, -1] = 1
    M[0, :] = 0
    M[-1, :] = 0
    M[0, 0] = 1
    M[-1, -1] = 1

    # Simulation des déplacements
    simulation = np.zeros((n_steps, n_nodes))  # Déplacements
    simulation[0, :] = perturbation_initiale(x, L, position_pincement, hauteur_max)


    # Initialiser la vitesse initiale (optionnel)
    simulation[1, :] = simulation[0, :]

    # Simulation temporelle avec méthode des différences centrales
    for time in range(1, n_steps - 1):
        A = np.linalg.inv(M) @ (-K @ simulation[time])
        simulation[time + 1, :] = 2 * simulation[time, :] - simulation[time - 1, :] + (dt**2) * A
        for space in range(1, n_elements):
            r_position = (r**2 + (L/n_nodes * space)**2)**(1/2) # position relative de où a lieu le mouvement
            index =int((time- r_position/c0))
            pression_in_time[time] += p0 * c0 / (4*math.pi) *(simulation[index  + 1, space ] - simulation[index, space])/dt/r_position #in pascal 

    return simulation, pression_in_time 


# Initialiser la perturbation
# Initialiser la perturbation avec une hauteur maximale de 0.01 m
hauteur_max = 0.0003  # Par exemple, 1 cm
#%% Initialisation
r = 1 #distance to where we are listening
p0 = 1.225 #ρ0​ : densité de l'air (∼1.225 kg/m3∼1.225kg/m3).
c0 = 343 #c0​ : vitesse du son dans l'air (∼343 m/s∼343m/s).
# Paramètres de la corde
L = 0.655              # Longueur de la corde (m)

position_pincement = 21  # La corde est pincée à un quart de sa longueur
rho = 1150*(0.00069/2)**2* np.pi       # Densité linéique de masse (kg/m)
T = 42.86                 # Tension de la corde (N)
c = (T / rho) ** 0.5   # Célérité de l'onde
print("Célérité de l'onde : c =", c)
# Condition CFL pour le pas de temps
t_end = 1         # Temps de simulation (s)
dt = 0.0000151
n_elements = 79    # Nombre d'éléments finis
dx = L / n_elements    # Longueur d'un élément
n_nodes = n_elements + 1
n_steps = int(t_end//dt)        # Nombre de pas temporels
start_time = time.time()
simulation, pression = FEMgui(L, rho, T,c, n_elements, dx, n_nodes, position_pincement, hauteur_max, t_end, dt, n_steps)
end_time = time.time()


elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time} seconds.")

# Enregistrer les données
sample_rate = int(1 / dt)
audio_data = np.int16((pression / np.max(np.abs(pression))) * 32767)
print(sample_rate)
file_name = f"el_finisimulation_L{L}_p0{p0}_c0{c0}_celerity{c}.wav"
with wave.open(file_name, 'w') as wf:
    wf.setnchannels(1)
    wf.setsampwidth(2)
    wf.setframerate(sample_rate)
    wf.writeframes(audio_data.tobytes())

# Affichage des résultats
plt.figure(1)
plt.plot(pression)
plt.ylabel('Pression (Pa)')
plt.xlabel('Temps (s)')
plt.title('Pression dans le temps à 1 m')
plt.show()


# Définir l'axe des positions (pour la simulation spatiale)
X = np.linspace(0, L, n_nodes)

# Tracé de la hauteur de la corde pour différents temps
temps = [int(k * n_steps / 30) for k in range(30)]  # Choisir 30 moments uniformément espacés
plt.figure(figsize=(15, 10))
for t in temps:
    plt.plot(X, simulation[t, :], label=f"t = {t * dt:.3f} s")
plt.ylabel('Hauteur (m)')
plt.xlabel('Position (m)')
plt.title('Hauteur de la corde en fonction de x pour différents temps')
plt.legend()
plt.show()
