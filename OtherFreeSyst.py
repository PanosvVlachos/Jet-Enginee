import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

plt.rcParams['figure.figsize'] = (12, 12)  # Use tuple () instead of set {}
plt.rcParams.update({'font.size': 18})     # Use .update(), not .updates()

a      = 1 # Diffusion
psi    = 0.5 # Distabilization term
L      = 2*np.pi
N_grid = 1000 #Number of discretization points
dx     = L/N_grid
x      = np.arange(-L/2, L/2, dx) #Define x domain

#Define descrete wavenumbers
kappa  = 2*np.pi*np.fft.fftfreq(N_grid, d= dx)

#Initial Conditions
u_0    = np.zeros_like(x)
u_0    = np.sin(x) #Ask AI for initial cond
u0hat  = np.fft.fft(u_0)
 #SciPy's odeint doesn't play well with complex number. 
 
u0hat_ri = np.concatenate((u0hat.real, u0hat.imag))

#Simulate in Fourier Frequency domain
dt     = 0.1
t      =np.arange(0, 10, dt)

def rhsEq(uhat_ri, t, kappa, a):
    uhat     = uhat_ri[:N_grid] + (1j) * uhat_ri[N_grid:]
    d_uhat    = -a**2 * (np.power(kappa,2) +  psi) * uhat 
    d_uhat_ri = np.concatenate((d_uhat.real, d_uhat.imag)).astype('float64') 
    return d_uhat_ri

uhat_ri   = odeint(rhsEq, u0hat_ri, t, args =(kappa, a))

uhat      = uhat_ri[:, :N_grid] + (1j) * uhat_ri[:, N_grid:]    
 
u         = np.zeros_like(uhat)

for k in range(len(t)):
    u[k, :] = np.fft.ifft(uhat[k, :])
    
u     =u.real

#Waterfall plot
fig = plt.figure()
ax  =fig.add_subplot(111, projection = '3d')
plt.set_cmap('jet_r')
u_plot = u[0:-1:10, : ]
for j in range(u_plot.shape[0]):
    ys = j*np.ones(u_plot.shape[1])
    ax.plot(x, ys, u_plot[j, :],color = cm.jet(j*20))


#Image plot
plt.figure()
plt.imshow(np.flipud(u), aspect=8)
plt.axis('off')
plt.set_cmap('jet_r')
plt.show()
#%%

plt.figure(figsize=(10, 6))

# Plot every Nth time slice to avoid clutter
for t in range(0, u.shape[0], 15):
    plt.plot(x, u[t, :], label=f'Time step {t}')

plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Classic Line Plot of u(x) at Different Time Steps')
plt.legend()
plt.grid(True)
plt.show()






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    