import numpy as np
import scipy.io.wavfile as wavfile

# Parameters
N = 90  # Number of points in theta
K = 180  # Number of points in phi
n0 = 40  # Point max in theta <= N
k0 = 100  # Point max in phi <= K
u0_i = 0.0001  # Amplitude at impact point at t0
Nt = 44100  # Number of points in time
Tau = 0.01  # Signal duration in seconds
R = 0.05  # Radius of the drum (m)
h = 0.001  # Thickness of the drum (m)
rho = 8500.0  # Density (kg/m^3)
sigma = 0.001  # Damping (kg/m^2/s)
E = 0  # Young's modulus
nu = 0.37  # Poisson coefficient
Gain = 100.0  # Amplification gain

# Filename
Lenom_Sce = 'Son_Timbre_Test.wav'

# Initial and boundary calculations
D = E * h**3 / 12 / (1 - nu**2)  # Flexural rigidity
delta_theta = np.pi / 2 / N  # Step in theta (radians)
delta_phi = 2 * np.pi / K  # Step in phi (radians)
delta_t = Tau / Nt  # Time step (seconds)

alpha = rho * h / delta_t**2 + sigma / (2 * delta_t)
gamma = 2 * rho * h / delta_t**2
zeta = -rho * h / delta_t**2 + sigma / (2 * delta_t)

theta = np.arange(1, N + 1) * delta_theta  # Antilatitude angles
phi = np.arange(1, K + 1) * delta_phi  # Longitude angles

# Initialize elongations
u1 = np.zeros((N, K))
u0 = np.zeros((N, K))
u_ = np.zeros((N, K))
Fs = Nt / Tau  # Sampling frequency
y_moy = np.zeros(Nt)  # Average level

# Initial shape of the drum
u0[n0, k0] = u0_i

# Main time loop
for t in range(1, Nt + 1):
    t_courant = t * delta_t
    print(t_courant)

    y_moy[t - 1] = Gain * np.mean(u0)

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

    for k in range(K):
        for n in range(1, N - 1):
            d_dtheta[n, k] = 1 / (2 * delta_theta) * (u0[n + 1, k] - u0[n - 1, k])
            d2_dtheta2[n, k] = 1 / delta_theta**2 * (u0[n + 1, k] - 2 * u0[n, k] + u0[n - 1, k])

    for k in range(1, K - 1):
        for n in range(N):
            d2_dphi2[n, k] = 1 / delta_phi**2 * (u0[n, k + 1] - 2 * u0[n, k] + u0[n, k - 1])

    # Special boundary cases
    for k in [0, K - 1]:
        for n in range(N):
            d2_dphi2[n, k] = 1 / delta_phi**2 * (u0[n, (k + 1) % K] - 2 * u0[n, k] + u0[n, (k - 1) % K])

    # Higher-order derivatives
    for k in range(K):
        for n in range(2, N - 2):
            d3_dtheta3[n, k] = 3 / (8 * delta_theta**3) * (u0[n + 2, k] - u0[n - 2, k] - 2 * (u0[n + 1, k] - u0[n - 1, k]))
            d4_dtheta4[n, k] = 1 / delta_theta**4 * (u0[n + 2, k] + u0[n - 2, k] - 4 * u0[n + 1, k] - 4 * u0[n - 1, k] + 6 * u0[n, k])

    for k in range(2, K - 2):
        for n in range(N):
            d4_dphi4[n, k] = 1 / delta_phi**4 * (
                u0[n, k + 2] + u0[n, k - 2] - 4 * u0[n, k + 1] - 4 * u0[n, k - 1] + 6 * u0[n, k])

    # Compute double Laplacian
    for k in range(K):
        for n in range(N):
            nabla4[n, k] = (
                d4_dtheta4[n, k]
                + 2 * np.cos(theta[n]) * d3_dtheta3[n, k]
                + ((np.cos(theta[n])**2 - 2) / (np.sin(theta[n])**2)) * d2_dtheta2[n, k]
                + (np.cos(theta[n]) / (np.sin(theta[n])**2)) * d_dtheta[n, k]
                + (1 / (np.sin(theta[n])**4)) * d4_dphi4[n, k]
                + 2 * ((np.cos(theta[n])**2 + 1) / (np.sin(theta[n])**4)) * d2_dphi2[n, k]
            )

    # Update elongations
    for k in range(K):
        for n in range(N):
            u1[n, k] = 1 / alpha * (-D / R**4 * nabla4[n, k] + gamma * u0[n, k] + zeta * u_[n, k])

    # Move to the next time step
    u_ = np.copy(u0)
    u0 = np.copy(u1)

# Write to WAV file
wavfile.write(Lenom_Sce, int(Fs), y_moy.astype(np.float32))
