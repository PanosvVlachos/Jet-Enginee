# Jet-Enginee
This Project is part of designing a controller for the Viscous Moore-Greitzer Model, a model that captures the dynamical behavior of a compression system.The model consists of a coupled system comprising one PDE, which describes the behavior of disturbances in the inlet region of compression systems, and two ODEs, which represent the coupling between these disturbances and the mean flow.

$$
\small
\frac{\partial}{\partial t} \begin{bmatrix}
    g\\
    \Phi\\
    \Psi
\end{bmatrix}= 
\begin{bmatrix}
    K^{-1} \left( \frac{\nu}{2} \frac{\partial^2}{\partial \theta^2} - \frac{1}{2}\frac{\partial}{\partial \theta} \right) & 0 & 0 \\
    0 & 0 & 0 \\
    0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
    g\\
    \Phi\\
    \Psi
\end{bmatrix}+
\begin{bmatrix}
    \alpha K^{-1} \left( \psi_c(\Phi+g)-\bar{\psi}_c(\Phi,g)\right) \\
    \frac{1}{l_c}\left(\bar{\psi}_c(\Phi,g)-\Psi\right) \\
    \frac{1}{4l_cB^2}\left(\Phi-\mu\sqrt{\Psi} \right) 
\end{bmatrix}  \color{black}.
\normalsize
$$

By linearizing the system around its equilibrium point, we observe that it becomes fully decoupled. To analyze the system dynamics, we apply the semigroup framework and use the spectral decomposition method—a classical yet powerful tool—to transform the PDE operator into an infinite set of ordinary differential equations (ODEs). This reformulation simplifies the control design process. A feedback control law is then developed using the Linear Quadratic Regulator (LQR) approach. To assess stability, a suitable Lyapunov functional is constructed handling the spillover effect, ensuring the system's robust performance over time. 


