import sys

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

target = r"""\begin{figure}[H]
    \centering
    \begin{subfigure}[b]{0.48\textwidth} \centering \includegraphics[width=\textwidth]{Images/fdtd_snapshots/swr_parameters/swr_empirical_error_2d.png} \caption{Measured empirical error (2D).} \end{subfigure}
    \hfill
    \begin{subfigure}[b]{0.48\textwidth} \centering \includegraphics[width=\textwidth]{Images/fdtd_snapshots/swr_parameters/swr_predicted_rho_2d.png} \caption{Predicted max convergence factor (2D).} \end{subfigure}
    \caption{Comparison of the empirically measured optimum error versus the theoretical continuous optimum over the parameter space $(q,r)$ for the 2D problem. The inclusion of the transverse $k_y$ spectrum successfully aligns the theoretical prediction with the numerical optimum.}
    \label{fig:oswr_empirical_vs_theoretical_2d}
\end{figure}"""

replacement = target + r"""

Furthermore, Figure~\ref{fig:oswr_convergence_comparison_2d} tracks the global $L_\infty$ error as a function of the Schwarz iteration count. It compares the convergence history of the unoptimized classical Robin transmission conditions ($q=1/c, r=0$) against the optimally tuned parameters $q^\ast$ and $r^\ast$. Utilizing the theoretically optimized transmission boundary operators dramatically accelerates the solver convergence, dropping the global error significantly faster and validating the efficacy of the 2D parameter optimization.

\begin{figure}[H]
    \centering
    \includegraphics[width=0.75\textwidth]{Images/fdtd_snapshots/swr_parameters/swr_convergence_comparison_2d.png}
    \caption{Comparison of the convergence rates of the SWR iteration for the 2D wave equation between unoptimized classical Robin parameters and the theoretically optimized configuration. The optimized parameters rapidly contract the iteration error towards the spatial discretization floor.}
    \label{fig:oswr_convergence_comparison_2d}
\end{figure}"""

tex = tex.replace(target, replacement)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

