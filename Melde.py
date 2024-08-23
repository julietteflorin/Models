"""
Created on Fri Jan 22 14:55:43 2021

@author: Juliette FLORIN
"""

# Importation des librairies utiles

import numpy as np
import matplotlib.pyplot as plt

def EqMelde(pulse, celerity, length, length_points, time_differential, time_points):
    simulation_tensor = np.zeros((length_points,time_points))
    x_differential = length / length_points
    
    # Conditions initiales
    for position in range(0, length_points):
        simulation_tensor[position, 0] = 0
        simulation_tensor[position, 1] = 0
    
    # Condition aux limites
    for time in range(1, time_points):
        simulation_tensor[0, time] = np.sin(pulse * time * time_differential) / 1000 # mise en mouvement transverse
        simulation_tensor[length_points - 1, time] = 0 # f=Fixe de l'autre côté
    
    # Calcul au cours du temps
    for time in range(1, time_points - 1):
        for position in range(1, length_points - 1):
            simulation_tensor[position,time + 1] = (-simulation_tensor[position, time - 1]
                                                    + 2 * simulation_tensor[position, time]
                                                    +(simulation_tensor[position - 1, time] 
                                                    + simulation_tensor[position + 1, time] 
                                                    - 2 * simulation_tensor[position, time]) * (celerity * time_differential / x_differential) ** 2
                                                    )
            
    return simulation_tensor
#%% Initialisation
pulse = 209.4
celerity = 10 # la célérité : 10m/s
length = 0.3
length_points = 100
time_differential =  10 ** (-5) 
time_points =  10 ** 5 # Pendant 1 seconde
#%% Simulation
simulation = EqMelde(pulse, celerity, length, length_points, time_differential, time_points) 


#%% Affichage tout sur une figure:

X = np.linspace(0, 0.3, 101)
 

plt.figure(4, (25, 25))
plt.xlim(0, 0.36) # Je vais plus loin, pour que les légendes ne gènent pas !
plt.ylim(-0.01, 0.01)
temps=[1000 * k for k in range(31)]
for t in temps:
    plt.plot(X, simulation[:, t], label = "tps : "+str(t // 10) + " ms.")
plt.ylabel('hauteur en m') # Légendes en x et y
plt.xlabel('x en m')
plt.title('hauteur de la corde en fonction de x pour différents temps')
plt.legend(loc = 7)
plt.show()




