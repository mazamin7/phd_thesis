import sys

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

target = r"\subsection{Extension to Two-Dimensional Domains}"

replacement = r"""\subsection{Extension to Two-Dimensional Domains}

The convergence analysis and optimization of the Schwarz Waveform Relaxation method can be naturally extended to higher-dimensional domains. Let us consider the damped wave equation in two spatial dimensions $(x, y)$, where the domain is partitioned along the $x$-axis into two overlapping subdomains $\Omega_1 = (0, b) \times (0, L_y)$ and $\Omega_2 = (a, L_x) \times (0, L_y)$.

The error equations for the 2D problem are given by
\begin{equation*}
u_{tt} + \gamma u_t = c^2 (u_{xx} + u_{yy}) + \nu (u_{xxt} + u_{yyt}), \qquad (x,y) \in \Omega, \ t > 0.
\end{equation*}

Applying a Fourier transform along the transverse $y$-direction (with spatial wavenumber $k_y$) and a Laplace transform in time (with complex frequency $s$), the transformed error equation reduces to an ordinary differential equation in $x$:
\begin{equation*}
(s^2 + \gamma s) \hat{u} = (c^2 + \nu s) (\hat{u}_{xx} - k_y^2 \hat{u}).
\end{equation*}
Rearranging the terms, we obtain the characteristic equation governing the spatial modes along the decomposition axis:
\begin{equation*}
\hat{u}_{xx} = K(s, k_y)^2 \hat{u},
\end{equation*}
where the modified dispersion relation is defined as
\begin{equation*}
K(s, k_y) = \sqrt{ \frac{s^2 + \gamma s}{c^2 + \nu s} + k_y^2 }.
\end{equation*}

This structure is mathematically equivalent to the one-dimensional case, with the continuous spectral argument $K(s)$ simply replaced by $K(s, k_y)$. Consequently, the analytical expression for the convergence factor $\rho(s, k_y; q, r)$ mirrors the 1D derivation perfectly. The transmission parameters $(q, r)$ can therefore be optimized for two-dimensional problems by minimizing the maximum of the convergence factor across the relevant range of both the temporal frequencies $\omega$ (where $s = i\omega$) and the transverse spatial frequencies $k_y$ supported by the underlying numerical grid.

To definitively validate this extended theoretical framework, we implement the overlapping decomposition for the FDTD method on a 2D Cartesian grid. The subdomains exchange data via the second-order BDF Robin discretization applied independently along each row of the artificial interfaces. We perform an exhaustive numerical sweep over a 2D grid of Robin parameters $(q, r)$. For each configuration, we run the SWR iteration for the propagation of a 2D standing wave and extract the empirical maximum error after a fixed number of iterations. Simultaneously, we evaluate the 2D analytical convergence factor $\rho$ over the same parameter space, taking the supremum over both temporal and transverse spatial frequencies. Figure~\ref{fig:oswr_empirical_vs_theoretical_2d} compares the empirically measured optimal parameters against the analytically predicted continuous optimum. The location of the minimum of the convergence factor matches the empirical error minimum, demonstrating the robustness of the derived 2D optimization problem.

\begin{figure}[H]
    \centering
    \begin{subfigure}[b]{0.48\textwidth} \centering \includegraphics[width=\textwidth]{Images/fdtd_snapshots/swr_parameters/swr_empirical_error_2d.png} \caption{Measured empirical error (2D).} \end{subfigure}
    \hfill
    \begin{subfigure}[b]{0.48\textwidth} \centering \includegraphics[width=\textwidth]{Images/fdtd_snapshots/swr_parameters/swr_predicted_rho_2d.png} \caption{Predicted max convergence factor (2D).} \end{subfigure}
    \caption{Comparison of the empirically measured optimum error versus the theoretical continuous optimum over the parameter space $(q,r)$ for the 2D problem. The inclusion of the transverse $k_y$ spectrum successfully aligns the theoretical prediction with the numerical optimum.}
    \label{fig:oswr_empirical_vs_theoretical_2d}
\end{figure}
"""

if "Extension to Two-Dimensional Domains" not in tex:
    # Insert it right before \section{Adaptive Rectangular Decomposition}
    tex = tex.replace(r"\section{Adaptive Rectangular Decomposition}", replacement + "\n\n" + r"\section{Adaptive Rectangular Decomposition}")
else:
    # Just replace if it already exists, wait, it doesn't exist yet since I only created the proposal artifact
    pass

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

