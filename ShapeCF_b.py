import numpy as np
# Define the shape function
def build_b_hat(a, k, psi, delta=0.5):
    N_full = len(k)
    unstable_k = []
    growth_rates = []

    # Identify unstable modes
    for i in range(N_full):
        ki = k[i]
        
        growth_rate = -a * ki**2 + psi
        if growth_rate > 0:
            unstable_k.append(ki)
            growth_rates.append(growth_rate)

    N_u = len(unstable_k)
    B = np.eye(N_u, dtype=complex)  # Identity matrix in unstable subspace
    u_vector = -np.array(growth_rates) - delta

    return B, u_vector, np.array(unstable_k)

def build_u_hat_full_dim(a, k, psi, delta=0.5):
    N = len(k)
    u_hat = np.zeros(N, dtype=complex)

    for i in range(N):
        ki = k[i]
        growth_rate = -a * ki**2 + psi
        if growth_rate > 0:
            u_hat[i] = -growth_rate - delta

    return u_hat