import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from mpl_toolkits.mplot3d import Axes3D

# Parameters
N = 128
x = np.linspace(0, 2*np.pi, N, endpoint=False)
dx = x[1] - x[0]
a = 0.4
b = 0.5
psi = 2
T = 2
dt = 0.0005
steps = int(T / dt)

# Initial condition
g0 = np.sin(x)

# Spectral multiplier
k = np.fft.fftfreq(N, d=dx) * 2 * np.pi
lam = -a * k**2 - 1j * b * k + psi

# RK4 step
def rk4_step(g_hat, dt, lam):
    k1 = lam * g_hat
    k2 = lam * (g_hat + 0.5 * dt * k1)
    k3 = lam * (g_hat + 0.5 * dt * k2)
    k4 = lam * (g_hat + dt * k3)
    return g_hat + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

# Time evolution
g_hat = fft(g0)
g_xt = [np.real(ifft(g_hat))]
for _ in range(steps):
    g_hat = rk4_step(g_hat, dt, lam)
    g_xt.append(np.real(ifft(g_hat)))
g_xt = np.array(g_xt)

# Create 3D surface plot
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# Time and space meshgrids
t_vals = np.linspace(0, T, steps + 1)
T_grid, X_grid = np.meshgrid(t_vals, x, indexing='ij')  # shapes match g_xt

# Plot the surface
surf = ax.plot_surface(X_grid, T_grid, g_xt, cmap='viridis', edgecolor='none')

# Labels and title
ax.set_xlabel("θ")
ax.set_ylabel("Time t")
ax.set_zlabel("g(θ, t)")
ax.set_title("3D Surface Plot of PDE Solution")

# Colorbar
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

plt.tight_layout()
plt.show()
