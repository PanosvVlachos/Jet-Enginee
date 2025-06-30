import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker


# Parameters
N_grid = 1000  # Fine grid
N_modes =16  # Number of Fourier modes to keep
x = np.linspace(0, 2*np.pi, N_grid, endpoint=False) #The spatial grid (1D array of length N)
n_unst = 4 #unstable states/ choose a suitable one(Both positive and negative frequencies)


def control_shape(n_unst, width=0.5, normalize=True):
    """
    Create `n_unst` localized shape functions over [0, 2π], sampled at N points.
    
    Parameters:
    - n_unst: Number of control points / shape functions
    - N_grid: Number of points to evaluate shape functions on (default 1000)
    - width: Width around the control point where the shape function is 1
    
    Returns:
    - b_shape: (n_unst, N) matrix of shape functions
    - x: The spatial grid (1D array of length N)
    - control_pnts: Locations of the control points (1D array of length n_unst)
    """
    control_pnts = np.linspace(0, 2 * np.pi, n_unst, endpoint=False)
    b_shape = np.zeros((n_unst, N_grid))
    for i, cp in enumerate(control_pnts):
        dist = np.abs((x - cp + np.pi) % (2 * np.pi) - np.pi)
        b_shape[i] = (dist < width / 2).astype(float)
        #if normalize:
         #   b_shape[i] /= np.sum(b_shape[i]) * (2 * np.pi / N_grid)  # Unit integral/ \int b^i\,d\theta=1
    return b_shape, control_pnts
   


#%%
b_shape, control_pnts = control_shape(n_unst)

for i in range(len(control_pnts)):
     plt.plot(x, b_shape[i], label=f"$b^{i}$")
plt.legend()
plt.title("Localized Shape Functions")
plt.xlabel(r"$\theta$")
plt.ylabel("Shape value $b^i$")
plt.grid(True)
# Customize x-ticks to show symbolic pi values
xticks = np.linspace(0, 2*np.pi, 5)  # Adjust number of ticks as needed
xtick_labels = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
plt.xticks(xticks, xtick_labels)

plt.tight_layout()

plt.show()
#%%

# Step 1: Generate shape functions
b_shape, control_pnts = control_shape(n_unst)

# Step 2: FFT of each shape function
b_shape_fft = np.fft.fft(b_shape, axis=1)

# Step 3: Zero out high frequencies, keeping only the first N_modes

k = np.fft.fftfreq(N_grid, d=2*np.pi/N_grid)  # Frequency bins

# Keep only the lowest N_modes frequencies (both positive and negative)
mask = np.zeros(N_grid, dtype=bool)
sorted_freqs = np.argsort(np.abs(k))  # Sort frequencies by magnitude
mask[sorted_freqs[:N_modes]] = True   # Keep the first N_modes lowest frequencies

# Apply the mask to each shape function's FFT
b_shape_fft_filtered = np.zeros_like(b_shape_fft, dtype=complex)
for i in range(n_unst):
    b_shape_fft_filtered[i] = b_shape_fft[i] * mask

# Step 4: Inverse FFT to get reconstructed shapes
b_shape_recon = np.fft.ifft(b_shape_fft_filtered, axis=1).real
#%%
# Step 5: Plot original and reconstructed shape functions (one figure per shape)
for i in range(n_unst):
    plt.figure(figsize=(8, 4))
    plt.plot(x, b_shape[i], label="Original", linestyle="--", alpha=0.5)
    plt.plot(x, b_shape_recon[i], label="Reconstructed", linewidth=2)
    plt.title(f"Shape Function {i} Reconstruction with {N_modes} Fourier Modes")
    plt.xlabel("θ")
    plt.ylabel("Function Value")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
#%%
# Step 5: Plot original and reconstructed shape functions
colors_1 = ['blue', 'orange', 'green', 'red']
colors_2 = ['purple', 'brown', 'pink', 'gray']
for i in range(n_unst):
    color1 = colors_1[i % len(colors_1)]
    color2 = colors_2[i % len(colors_2)]

    # Original shape
    plt.plot(x, b_shape[i], color=color1)

    # Reconstructed shape
    plt.plot(x, b_shape_recon[i], color=color2, marker='o', linestyle='dashed', linewidth=0.5, markersize=0.12)
plt.title(f"Reconstruction of Shape Functions with {N_modes} Fourier Modes")
plt.xlabel(r"$\theta$")
plt.ylabel("Values of $b^i$")    

# Customize x-ticks to show symbolic pi values
xticks = np.linspace(0, 2*np.pi, 5)  # Adjust number of ticks as needed
xtick_labels = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
plt.xticks(xticks, xtick_labels)

plt.tight_layout()
plt.legend()
plt.grid(True)
plt.show()    

#%%
B = np.real(b_shape_fft[:, 1:n_unst+1]).T
print(B)
#%%
B_stabe = np.real(b_shape_fft[:,n_unst:N_modes]).T
print(B_stabe)
#%%
nu=0.5
C=1
# Initialize Xi matrix
Xi = 0.0

# Sum over strictly n > N (positive frequencies only)
for i in range(N_grid):
    if k[i] > n_unst:  # Strictly n > N condition
        denominator = (nu/2) * k[i]**2 - C
        bn = b_shape_fft[:, i]  # FFT coefficients at this mode
        Xi += np.dot(bn.conj(), bn).real / denominator  # Dot product gives scalar
#%%