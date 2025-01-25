import numpy as np
import matplotlib.pyplot as plt
import time
import wave
import math
def calculate_string_vibration(N, n0, y0, K, Tau, L, T, sigma, mu, delta_x, delta_t):
    """
    Calculate the vibration of a string and return the mean level of vibration.

    Parameters:
        N (int): Number of points in space
        n0 (int): Index of the maximum point, must be <= N
        y0 (float): Maximum amplitude
        K (int): Number of points in time
        Tau (float): Duration of the signal in seconds
        L (float): Length of the string in meters
        T (float): Tension of the string in Newtons
        sigma (float): Damping in kg/m/s
        mu (float): Linear mass in kg/m
        Gain (float): Gain amplification on the emitted level

    Returns:
        y (numpy.ndarray): Displacement of the string over time
        pression_in_time (numpy.ndarray): Pressure over time at a fixed position
    """
    # Initialize parameters
    p0 = 1.225 # ρ0: density of air (~1.225 kg/m^3).
    c0 = 343 # c0: speed of sound in air (~343 m/s)
    gamma = T * delta_t**2 / delta_x**2
    alpha = mu + sigma / 2 * delta_t
    theta = -mu + sigma / 2 * delta_t
    r = 1
    y = np.zeros((N, K))
    pression_in_time = np.zeros(K)
    
    # Shape of the string and initial velocities
    for n in range(0, n0):
        y[n, 0] = y0 *n / n0
    for n in range(n0, N):
        y[n, 0] = y0 * (N - n) / (N - n0)

    for time in range(1, K):
        y[0, time] = 0 
        y[N - 1, time] = 0
    for position in range(0, N):
        y[position, 1] = y[position, 0]
    # Current values calculation
    for k in range(1, K-1):
        for n in range(1, N - 1):
            y[n, k + 1] = (gamma * (y[n + 1, k] + y[n - 1, k]) +
                           2 * (mu - gamma) * y[n, k] +
                           theta * y[n, k - 1]) / alpha
            r_position = (r**2 + (L/N * n)**2)**(1/2) # position relative de où a lieu le mouvement
            index = int((k - r_position / c0))
            pression_in_time[k] += p0 * c0 / (4 * math.pi) * ( y[n, index + 1] - y[n, index] ) / delta_t / r_position

    return y, pression_in_time
    

# Parameters
N = 80
n0 = 22
y0 = 0.0003
diff_t = 1*10**(-5)
Tau = 7
K = int(Tau//diff_t)

L = 0.655
T = 42.86 # kg into Newton
sigma = 0
mu = 1150*(0.00069/2)**2* np.pi #0.000582
print(mu)
delta_x = L / N
#%% Simulation
start_time = time.time()
simulation, pression = calculate_string_vibration(N, n0, y0, K, Tau, L, T, sigma, mu, delta_x, diff_t)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time} seconds.")

#%% Save the data
file_name = (
    f"comparing_simulation_N{N}_n0{n0}_y0{y0}_Tau{Tau}_L{L}_T{T}_sigma{sigma}_mu{mu}.wav"
)
file_name = file_name.replace(" ", "_")

sample_rate = 1 / diff_t
num_samples = len(pression)
amplitude = 32767

pression = pression / np.max(np.abs(pression))
audio_data = np.int16(pression * amplitude)

with wave.open(file_name, 'w') as wf:
    wf.setnchannels(1)
    wf.setsampwidth(2)
    wf.setframerate(sample_rate)
    wf.writeframes(audio_data.tobytes())

#%% Plotting
X = np.linspace(0, L, N)

plt.figure(figsize=(12, 8))
times_to_plot = [k for k in range(31)]
delta_t = K//31
for t in times_to_plot:
    plt.plot(X, simulation[:, t*delta_t], label=f"Time: {(t * delta_t)*Tau/K:.3f} s")
plt.ylabel('Height (m)')
plt.xlabel('Position (m)')
plt.title('String height vs position for different times')
plt.legend(loc='best')
plt.show()

plt.figure(figsize=(12, 8))
plt.plot(pression)
plt.ylabel('Pressure (Pa)')
plt.xlabel('Time (s)')
plt.title('Pressure over time at 1 m')
plt.show()
