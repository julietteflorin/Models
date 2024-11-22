"""
Created on Fri Jan 22 14:55:43 2021

@author: Juliette FLORIN
"""

# Importation des librairies utiles

import numpy as np
import matplotlib.pyplot as plt
import math
import time
import wave
def Eqguitard(r, p0, c0, U_0, celerity, length, length_points, time_differential, time_points, position_tire):
    pression_in_time = np.zeros(time_points)
    
    simulation_tensor = np.zeros((length_points,time_points))
    #pression_accoustique = p0 * c0 * vitesse transversale instantanée de la corde * 1/ distance cherché 
    x_differential = length / length_points
    
    # Conditions initiales
    for position in range(0, position_tire):
        simulation_tensor[position, 0] = U_0 * position/position_tire
    for position in range(position_tire, length_points):    
        simulation_tensor[position, 0] = U_0 * (length - position * (length/length_points) /(length- position_tire** (length/length_points)))
    
    # Calcul au cours du temps----
    for time in range(1, time_points - 1):
        for position in range(1, length_points - 1):
            simulation_tensor[position,time + 1] = (-simulation_tensor[position, time - 1]
                                                    + 2 * simulation_tensor[position, time]
                                                    +(simulation_tensor[position - 1, time] 
                                                    + simulation_tensor[position + 1, time] 
                                                    - 2 * simulation_tensor[position, time]) * (celerity * time_differential / x_differential) ** 2
                                                    )
            r_position = (r**2 + (length/length_points * position)**2)**(1/2) # position relative de où a lieu le mouvement
            index =int((time- r_position/c0))
            pression_in_time[time] += p0 * c0 / (4*math.pi) *(simulation_tensor[position, index  + 1] - simulation_tensor[position, index])/time_differential/r_position 
    return simulation_tensor, pression_in_time
#%% Initialisation
r = 1 #distance to where we are listening
p0 = 1.225 #ρ0​ : densité de l'air (∼1.225 kg/m3∼1.225kg/m3).
c0 = 343 #c0​ : vitesse du son dans l'air (∼343 m/s∼343m/s).
celerity = 2600 # la célérité : 10m/s
length = 0.5
length_points = 200
time_differential =  10 ** (-7)#-6 error 
time_points =  10 ** 5 #1 ms
U_0 = 0.0003
position_tire = length_points//3
#%% Simulation
start_time = time.time()
simulation, pression = Eqguitard(r, p0, c0, U_0, celerity, length, length_points, time_differential, time_points, position_tire) 
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time} seconds.")

#%% Enregistrer les données
file_name = (
    f"simulation_r{r}_p0{p0}_c0{c0}_celerity{celerity}_"
    f"length{length}_points{length_points}_dt{time_differential}_"
    f"timepoints{time_points}_U0{U_0}.wav"
)
sample_rate = len(pression)/(len(pression)*time_differential)
num_samples = len(pression)
amplitude = 32767

pression = pression / np.max(np.abs(pression))
audio_data = np.int16(pression * amplitude)

with wave.open(file_name, 'w') as wf:
    wf.setnchannels(1)
    wf.setsampwidth(2)
    wf.setframerate(sample_rate)
    wf.writeframes(audio_data.tobytes())
#%% Upload les données
#%% Affichage tout sur une figure:

X = np.linspace(0, length, length_points)

plt.figure(4, (25, 25))
temps=[1000 * k for k in range(31)]
for t in temps:
    plt.plot(X, simulation[:, t], label = "tps : "+str(t // 10) + " ms.")
plt.ylabel('hauteur en m') # Légendes en x et y
plt.xlabel('x en m')
plt.title('hauteur de la corde en fonction de x pour différents temps')
plt.legend(loc = 7)
plt.show()

plt.figure(4, (25, 25))
plt.plot(pression)
plt.ylabel('hauteur en m') # Légendes en x et y
plt.xlabel('x en m')
plt.title('pression dans le temps à 1 m')
plt.show()





