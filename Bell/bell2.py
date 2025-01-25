import numpy as np
import scipy.io.wavfile as wavfile
import matplotlib.pyplot as plt
import pickle 
import time

# Barre de progression pour suivre l'avancée de la simulation
def progress_bar(progress,totalIterations,start_time,bar_length=30):
    fraction = progress / (totalIterations - 1)
    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '

    elapsed_time = time.time() - start_time
    fraction = progress / (totalIterations - 1)
    remaining_time = elapsed_time/fraction - elapsed_time if fraction !=0 else 0

    hours, remainder = divmod(remaining_time, 3600)
    minutes, seconds = divmod(remainder, 60)

    p = ""
    if remaining_time != 0:
        p += f" {int(hours)} h" if int(hours) > 0 else ""
        p += f" {int(minutes)} m" if int(minutes) > 0 else ""
        p += f" {int(seconds)} s" if int(seconds) > 0 else ""
        p += " remaining"

    print('\r'+f'Progress: [{arrow}{padding}] {int(fraction*100)}%{p}'+' '*10, end="")

def simulation(R, u0_i):
    # Parameters
    N = 80  # Number of points in theta
    K = 160  # Number of points in phi
    n0 = 40  # Point max in theta <= N
    k0 = 100  # Point max in phi <= K
    Nt = 44100  # Number of points in time
    Tau = 0.01  # Signal duration in seconds
    h = 0.001  # Thickness of the drum (m)
    rho = 8500.0  # Density (kg/m^3)
    sigma = 0.001  # Damping (kg/m^2/s)
    E = 0  # Young's modulus
    nu = 0.37  # Poisson coefficient
    Gain = 100.0  # Amplification gain

    # Filename
    

    # Initial and boundary calculations
    D = E * h**3 / 12 / (1 - nu**2)  # Flexural rigidity
    delta_theta = np.pi / 2 / N  # Step in theta (radians)
    delta_phi = 2 * np.pi / K  # Step in phi (radians)
    delta_t = Tau / Nt  # Time step (seconds)



    theta = np.arange(1, N + 1) * delta_theta  # Antilatitude angles

    # Initialize elongations
    u1 = np.zeros((N, K))
    u0 = np.zeros((N, K))
    u_ = np.zeros((N, K, Nt), dtype=np.float32)
    try:
        u_ = np.zeros((N, K, Nt), dtype=np.float64)
    except Exception as e:
        print(e)
        print("Using dtype=np.float32 instead")

    Fs = Nt / Tau  # Sampling frequency
    y_moy = np.zeros(Nt)  # Average level

    # Initial shape of the drum
    u0[n0, k0] = u0_i
    u_[:,:,0]= u0

    # Derivatives
    d_dtheta = np.zeros((N, K))
    d2_dtheta2 = np.zeros((N, K))
    d2_dphi2 = np.zeros((N, K))
    d3_dtheta3 = np.zeros((N, K))
    d4_dtheta4 = np.zeros((N, K))
    d4_dphi4 = np.zeros((N, K))
    d3_dphi_dtheta2 = np.zeros((N, K))
    d3_dtheta_dphi2 = np.zeros((N, K))
    d4_dtheta2_dphi2 = np.zeros((N, K))
    nabla4 = np.zeros((N, K))

    # Main time loop
    start_time = time.time()
    for t in range(1, Nt-1):
        progress_bar(t,Nt,start_time)
        u0 = u_[:,:,t]
        y_moy[t - 1] = Gain * np.mean(u0)

        d_dtheta[1:N-1, :] = (u0[2:N, :] - u0[0:N-2, :]) / (2 * delta_theta)
        d2_dtheta2[1:N-1, :] = (u0[2:N, :] - 2 * u0[1:N-1, :] + u0[0:N-2, :]) / delta_theta**2

        d2_dphi2[:, 1:K-1] = (u0[:, 2:K] - 2 * u0[:, 1:K-1] + u0[:, 0:K-2]) / delta_phi**2

        # Higher-order derivatives
        d3_dtheta3[2:N-2, :] = (3 / (8 * delta_theta**3)) * (
            u0[4:N, :] - u0[0:N-4, :] - 2 * (u0[3:N-1, :] - u0[1:N-3, :])
        )
        d4_dtheta4[2:N-2, :] = (1 / delta_theta**4) * (
            u0[4:N, :] + u0[0:N-4, :] - 4 * u0[3:N-1, :] - 4 * u0[1:N-3, :] + 6 * u0[2:N-2, :]
        )
                
        d4_dphi4[:, 2:K-2] = (1 / delta_phi**4) * (
            u0[:, 4:K] + u0[:, 0:K-4] - 4 * u0[:, 3:K-1] - 4 * u0[:, 1:K-3] + 6 * u0[:, 2:K-2]
        )

        d3_dphi_dtheta2[1:N-1, 1:K-1] = (1 / (2 * delta_phi * delta_theta**2)) * (
            u0[2:N, 2:K] - 2 * u0[1:N-1, 2:K] + u0[0:N-2, 2:K] +
            u0[2:N, 1:K-1] + 2 * u0[1:N-1, 1:K-1] - u0[0:N-2, 1:K-1]
        )

        d3_dtheta_dphi2[1:N-1, 1:K-1] = (1 / (2 * delta_theta * delta_phi**2)) * (
            u0[2:N, 2:K] - 2 * u0[2:N, 1:K-1] + u0[2:N, 0:K-2] -
            u0[0:N-2, 2:K] + 2 * u0[0:N-2, 1:K-1] - u0[0:N-2, 0:K-2]
        )

        d4_dtheta2_dphi2[1:N-1, 1:K-1] = (1 / delta_theta**2) * (
            d2_dphi2[2:N, 1:K-1] - 2 * d2_dphi2[1:N-1, 1:K-1] + d2_dphi2[0:N-2, 1:K-1]
        )

        # Compute double Laplacian

        tan_theta = np.tan(theta)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        nabla4 = (
            d4_dtheta4 + 
            2 * (1 / tan_theta[:, None]) * d3_dtheta3 + 
            ((cos_theta**2 - 2) / (sin_theta**2))[:, None] * d2_dtheta2 + 
            (cos_theta / (sin_theta**2))[:, None] * d_dtheta + 
            (1 / (sin_theta**4))[:, None] * d4_dphi4 + 
            2 * ((cos_theta**2 + 1) / (sin_theta**4))[:, None] * d2_dphi2 + 
            2 * d4_dtheta2_dphi2 - 
            2 * (1 / tan_theta[:, None]) * d3_dphi_dtheta2
        )
        

        alpha = rho * h / delta_t**2 + sigma / (2 * delta_t)
        gamma = 2 * rho * h / delta_t**2
        zeta = -rho * h / delta_t**2 + sigma / (2 * delta_t)

        # Update elongations

        # Point fixe en haut, il reste à 0
        # Calcul pour tous les autres points sauf en pi/2
        u_[1:, :, t+1] = (1 / alpha) * (-D / R**4 * nabla4[1:, :] + gamma * u_[1:, :,t] + zeta * u_[1:, :, t-1])
        # Condition de liberté en pi/2
        u_[N-1, :, t + 1] = u_[N-2, :,t + 1]

    print()
    return u_, y_moy, Fs

# Pour l'affichage on crée une matrice
def create_display_tensor(simulation_tensor,n_display_points,radius):
    # On se place dans les coordonnées cartésiennes
    def to_cartesian_vector(amplitude,phi,theta):
        cartesian_array = [
            (radius + amplitude) * np.sin(theta) * np.cos(phi),
            (radius + amplitude) * np.sin(theta) * np.sin(phi),
            (radius + amplitude) * np.cos(theta)]
        return cartesian_array
    display_tensor = []
    for i in range(n_display_points): # par rapport au nombre de points qu'on veut afficher 
        time = i*int(np.shape(simulation_tensor)[0])//n_display_points
        positions = []
        for horizontal_angle in range(np.shape(simulation_tensor)[1]):
            phi = 2.*np.pi * horizontal_angle/np.shape(simulation_tensor)[1]
            for vertical_angle in range(np.shape(simulation_tensor)[0]):
                theta = np.pi/2. * vertical_angle/np.shape(simulation_tensor)[2]
                positions+=[to_cartesian_vector(simulation_tensor[time,horizontal_angle,vertical_angle],phi,theta)]
        display_tensor+=[positions]
    return display_tensor
# On a une matrice avec [[x,y,z],[x,y,z],...,[x,y,z]]
# On transpose cette matrice pour avoir [[x,...,x],[y,...,y],[z,...,z]]
def transpose(display_tensor):
    for i,time_point in enumerate(display_tensor):
        x,y,z = [],[],[]
        for position_vector in time_point:
            x += [position_vector[0]]
            y += [position_vector[1]]
            z += [position_vector[2]]
        display_tensor[i] = [x,y,z]
    return display_tensor
# Simulations pour différentes fréquences

Lenom_Sce = 'Son_Timbre_Test.wav'

n_display_points = 300
R = 0.05  # Radius of the drum (m)
u0_i = 0.0001  # Amplitude at impact point at t0
simulation_tensor, y_moy, Fs = simulation(R, u0_i)
# Write to WAV file
wavfile.write(Lenom_Sce, int(Fs), y_moy.astype(np.float32))
#%%  Mettre en coordonnées cartésiennes
display_tensor = create_display_tensor(simulation_tensor, n_display_points, R)
#%% Transposer la matrice coordonées cartésiens 
display_tensor_T = transpose(display_tensor)
#%% Enregistrer la simulation
my_file = 'simulation_freq_35449.9,5,40.txt'
with open(my_file, 'wb') as f:
    pickle.dump(simulation_tensor, f)

#%% Affichage
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

frames = len(display_tensor_T)
skip = 1
interval = 500 #Temps d'attente entre chaque cadre image

fig = plt.figure()
ax = Axes3D(fig)
def animate(i):
    ax.clear()
    i = skip*i # Intervalle de temps auquel on selectionne la partie à afficher
    x,y,z= (display_tensor_T[i][0],
            display_tensor_T[i][1],
            np.linalg.norm(np.array(display_tensor_T[i+1]) - np.array(display_tensor_T[i]),axis=0))
    # Afficher la sonnette
    #x,y,z= display_tensor_T[i][0], display_tensor_T[i][1],display_tensor_T[i][2] 
    ax.plot_trisurf(x,y,z,cmap=plt.cm.viridis, linewidth=0.2)
    ## On zoome en fonction de l'amplitude de la vitesse pour mieux voir
    ax.set_zlim([0,0.0056])                                                   
ani = animation.FuncAnimation(fig, 
                              animate, 
                              frames=int(frames), 
                              interval=interval, 
                              repeat=True)
plt.show()
