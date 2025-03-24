import numpy as np
import matplotlib.pyplot as plt
import time
import wave
import math
from scipy.linalg import solve
from numba import njit, prange

def perturbation_initiale(x, L, position_pincement, hauteur_max):
    perturbation = np.zeros_like(x)
    mid_point = int(position_pincement * len(x))  # Convertir en index
    for i in range(len(x)):
        if i <= mid_point:  # Pente ascendante
            perturbation[i] = hauteur_max * (i / mid_point)
        else:  # Pente descendante
            perturbation[i] = hauteur_max * ((len(x) - i) / (len(x) - mid_point))
    return perturbation

@njit(parallel=True)
def compute_pressure(simulation, r_positions, indices, dt, p0, c0, n_steps):
    pression_in_time = np.zeros(n_steps)
    for time in prange(1, n_steps - 1):
        for space in prange(1, len(r_positions)):
            index = indices[time, space]
            if 0 <= index < n_steps - 1:
                pression_in_time[time] += p0 * c0 / (4 * math.pi * r_positions[space]) * (
                    (simulation[index + 1, space] - simulation[index, space]) / dt
                )
    return pression_in_time

def FEMgui(L, rho, T, c, n_elements, dx, n_nodes, position_pincement, hauteur_max, t_end, dt, n_steps):
    x = np.linspace(0, L, n_nodes)
    dt_max = dx / c
    dt = min(dt, dt_max)
    print(f"dt utilisé : {dt}, dt_max (CFL) : {dt_max}")

    K_local = (T / dx) * np.array([[1, -1], [-1, 1]])
    M_local = (rho * dx / 6) * np.array([[2, 1], [1, 2]])
    K = np.zeros((n_nodes, n_nodes))
    M = np.zeros((n_nodes, n_nodes))
    for i in range(n_elements):
        dofs = [i, i + 1]
        K[np.ix_(dofs, dofs)] += K_local
        M[np.ix_(dofs, dofs)] += M_local

    K[0, :] = K[-1, :] = 0
    K[0, 0] = K[-1, -1] = 1
    M[0, :] = M[-1, :] = 0
    M[0, 0] = M[-1, -1] = 1

    simulation = np.zeros((n_steps, n_nodes))
    simulation[0, :] = perturbation_initiale(x, L, position_pincement, hauteur_max)
    simulation[1, :] = simulation[0, :]

    r_positions = np.sqrt(r**2 + (L / n_nodes * np.arange(1, n_elements))**2)
    indices = (np.arange(n_steps)[:, None] - r_positions / c0 / dt).astype(int)

    for time in range(1, n_steps - 1):
        A = solve(M, -K @ simulation[time])
        simulation[time + 1, :] = 2 * simulation[time, :] - simulation[time - 1, :] + (dt**2) * A

    pression_in_time = compute_pressure(simulation, r_positions, indices, dt, p0, c0, n_steps)
    return simulation, pression_in_time

position = 1
# Initialiser la perturbation
position_pincement = 0.25  # La corde est pincée à un quart de sa longueur
# Initialiser la perturbation avec une hauteur maximale de 0.01 m
hauteur_max = 0.0003  # Par exemple, 1 cm
#%% Initialisation
r = 1 #distance to where we are listening
p0 = 1.225 #ρ0​ : densité de l'air (∼1.225 kg/m3∼1.225kg/m3).
c0 = 343 #c0​ : vitesse du son dans l'air (∼343 m/s∼343m/s).
# Paramètres de la corde
L = 0.65               # Longueur de la corde (m)
rho = 0.000582         # Densité linéique de masse (kg/m)
T = 60                 # Tension de la corde (N)
c = (T / rho) ** 0.5   # Célérité de l'onde
print("Célérité de l'onde : c =", c)
# Condition CFL pour le pas de temps
t_end = 1         # Temps de simulation (s)
dt = 1*10**(-5)
n_elements = 80       # Nombre d'éléments finis
dx = L / n_elements    # Longueur d'un élément
n_nodes = n_elements + 1
n_steps = int(t_end//dt)        # Nombre de pas temporels
start_time = time.time()
simulation, pression = FEMgui(L, rho, T,c, n_elements, dx, n_nodes, position_pincement, hauteur_max, t_end, dt, n_steps)
end_time = time.time()


elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time} seconds.")


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

# Enregistrer les données
sample_rate = int(1 / dt)
audio_data = np.int16((pression / np.max(np.abs(pression))) * 32767)

file_name = f"simulation_r{r}_p0{p0}_c0{c0}_celerity{c}.wav"
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