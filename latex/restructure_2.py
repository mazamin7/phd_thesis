import re

with open("Thesis_restructured.tex", "r", encoding="utf-8") as f:
    tex = f.read()

target3 = r"where the modal amplitudes now include the effect of telegrapher and viscoelastic damping. This benchmark assesses the accuracy of the schemes for smooth damped standing waves. The FDTD method exhibits the expected second-order convergence of the discretization error. Since the exact solution is completely represented by the retained cosine modes, the PSTD method reproduces the solution up to machine precision."
repl3 = "where the modal amplitudes now include the effect of telegrapher and viscoelastic damping."
tex = tex.replace(target3, repl3)

target4 = r"The reference solution is computed from the cosine expansion of the initial condition, with each modal coefficient propagated analytically according to the damped modal dynamics. As in the undamped case, the modal expansion is truncated after a sufficiently large number of modes. This benchmark investigates the influence of damping on non-smooth initial data. The FDTD method is expected to retain second-order convergence of the truncation error. However, the limited spatial regularity of the initial profile introduces an initial error that propagates and cannot be completely distinguished from the discretization error, resulting in an observed convergence rate slightly below second order. For the PSTD method, the discretization error continues to converge algebraically."
repl4 = r"The reference solution is computed from the cosine expansion of the initial condition, with each modal coefficient propagated analytically according to the damped modal dynamics. As in the undamped case, the modal expansion is truncated after a sufficiently large number of modes."
tex = tex.replace(target4, repl4)

target5 = r"""\end{equation*}
This experiment models the propagation of a highly localized wave packet, introducing more spatial content compared to the previous experiments. Because the analytical solution consists of travelling waves, it provides a very rigorous test for the numerical dispersion of the FDTD method. Since the initial condition does not perfectly vanish at the boundaries, its true cosine spectrum does not decay exponentially, limiting the spatial accuracy of the PSTD method to algebraic convergence."""
repl5 = r"\end{equation*}"
tex = tex.replace(target5, repl5)

target6 = r"""because evaluating its cosine expansion would converge too slowly and introduce severe Gibbs phenomena. The final experiment evaluates the methods' capabilities in handling non-smooth propagating waves, emphasizing the effects of numerical dispersion on discontinuous profiles. Similar to the previous experiment, the PSTD method's convergence is constrained by the algebraic decay of the coefficients."""
repl6 = "because evaluating its cosine expansion would converge too slowly and introduce severe Gibbs phenomena."
tex = tex.replace(target6, repl6)

fdtd_eq = r"""\section{Finite-Difference Time-Domain Method}
In this section, we derive and analyze an FDTD scheme for the damped wave equation, following the approach described in \cite{bilbao_numerical_2009}, Chapter 6, \S6.2. For simplicity, we derive the method in one spatial dimension. The extension to higher dimensions is obtained by replacing the second spatial derivative by the corresponding discrete Laplacian \cite[Chapter 11]{bilbao_numerical_2009}.

The finite-difference time-domain (FDTD) method approximates the wave equation by replacing the continuous derivatives with finite difference operators on a discrete grid. For the damped wave equation, a standard explicit centered-difference scheme is given by:
\begin{equation}
\label{eq:fdtd_scheme}
\frac{u_i^{n+1} - 2u_i^n + u_i^{n-1}}{\Delta t^2}
- c^2 \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2}
+ \gamma \frac{u_i^{n+1} - u_i^{n-1}}{2\Delta t}
- \nu \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n - (u_{i+1}^{n-1} - 2u_i^{n-1} + u_{i-1}^{n-1})}{\Delta t \Delta x^2} = 0,
\end{equation}
where $u_i^n \approx u(i\Delta x, n\Delta t)$. This leads to an explicit update equation for $u_i^{n+1}$ in terms of the fields at time steps $n$ and $n-1$.
"""

tex = tex.replace(r"""\section{Finite-Difference Time-Domain Method}
In this section, we derive and analyze an FDTD scheme for the damped wave equation, following the approach described in \cite{bilbao_numerical_2009}, Chapter 6, \S6.2. For simplicity, we derive the method in one spatial dimension. The extension to higher dimensions is obtained by replacing the second spatial derivative by the corresponding discrete Laplacian \cite[Chapter 11]{bilbao_numerical_2009}.""", fdtd_eq)

pstd_eq = r"""\section{Pseudo-Spectral Time-Domain Method}
While the finite difference method evaluates spatial derivatives using local approximations on the grid, the Pseudo-Spectral Time-Domain (PSTD) method employs global orthogonal basis functions to evaluate spatial derivatives with spectral accuracy.

In PSTD, the spatial Laplacian is evaluated exactly up to the Nyquist limit using the Discrete Cosine Transform (DCT) to implicitly impose the homogeneous Neumann boundary conditions:
\begin{equation}
    \Delta u \approx \mathcal{C}^{-1} \left[ -\lambda_m \mathcal{C}[u] \right],
\end{equation}
where $\mathcal{C}$ denotes the DCT, and $\lambda_m = (m\pi/L)^2$ are the exact eigenvalues of the Laplacian operator on the discrete grid. The temporal integration is performed analytically by independently evolving each spectral mode using the exact state-transition matrix derived in Section~\ref{sec:modal_expansion}.
"""
tex = tex.replace(r"""\section{Pseudo-Spectral Time-Domain Method}
While the finite difference method evaluates spatial derivatives using local approximations on the grid""", pstd_eq)

with open("Thesis_restructured_2.tex", "w", encoding="utf-8") as f:
    f.write(tex)

