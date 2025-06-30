import numpy as np

# Define Gaussian shape centered at pi/2 and 3pi/2
def shape_fun(x, sigma):
    b_x = np.exp(-((x - np.pi/2)**2) / (2 * sigma**2)) + np.exp(-((x - 3*np.pi/2)**2) / (2 * sigma**2))

    #Normalize (optional)
    b_x /= np.max(b_x)
    return b_x
#%%
# Find unstable modes: lambda_k > 0, but next mode is stable
def number_of_unstable_modes(eigvals, k):
    unstable_modes = np.where(np.real(eigvals) > 0)[0]
    unstable_modes = unstable_modes[k[unstable_modes] != 0]
    unstable_k = k[unstable_modes]
    
    print("Then, the # of the unstable modes is:", len(unstable_k))
    print("Unstable modes:",sorted(unstable_k))
    return unstable_k
