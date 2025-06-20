# Jet-Enginee
This Project is part of designing a controller for the Viscous Moore-Greitzer Model, a model that captures the dynamical behavior of a compression system.The model consists of a coupled system comprising one PDE, which describes the behavior of disturbances in the inlet region of compression systems, and two ODEs, which represent the coupling between these disturbances and the mean flow.
<script type="text/javascript" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>
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
