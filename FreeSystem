import numpy as np

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from Engineering import number_of_unstable_modes

# Parameters
N = 128
x = np.linspace(0, 2*np.pi, N, endpoint=False)
dx = x[1] - x[0]
a = 0.7
b = 0.5
psi = 5
T = 5
dt = 0.0001
steps = int(T / dt)

# Initial condition
g0 =np.sin(x)

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

#Unstable modes
number = number_of_unstable_modes(lam, k)

#Control takes place

# Time indices to plot
snapshots = [0, steps // 3, 2 * steps // 3, steps]
times = [i * dt for i in snapshots]

# Plot snapshots
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True)
axes = axes.flatten()

for ax, idx, t in zip(axes, snapshots, times):
    ax.plot(x, g_xt[idx])
    ax.set_title(f"t = {t:.3f}")
    ax.set_xlabel("θ")
    ax.set_ylabel("g(θ, t)")
    ax.grid(True)

fig.suptitle("PDE Solution: ∂g/∂t = a gₓₓ - b gₓ + ψ'_(Φ_e) g", fontsize=14)
plt.tight_layout()
plt.show()
