import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

target = r"""\section{Pseudo-Spectral Time-Domain Method}
The modal decomposition presented above naturally suggests a pseudo-spectral spatial discretization. We follow the approach presented in \cite{raghuvanshi_efficient_2009}: as in the FDTD method, the PDE is enforced directly on a cell-centered Neumann grid, but rather than approximating the Laplace operator by finite differences, its action is evaluated spectrally by expanding the solution in the DCT basis, where the discrete Laplacian is diagonal. Finally, the evolution problem of the resulting uncoupled modal equations is solved using an exact time integrator."""

repl = r"""\section{Pseudo-Spectral Time-Domain Method}
The modal decomposition presented above naturally suggests a pseudo-spectral spatial discretization. We follow the approach presented in \cite{raghuvanshi_efficient_2009}: as in the FDTD method, the PDE is enforced directly on a cell-centered Neumann grid, but rather than approximating the Laplace operator by finite differences, its action is evaluated spectrally by expanding the solution in the DCT basis, where the discrete Laplacian is diagonal.

In PSTD, the spatial Laplacian is evaluated exactly up to the Nyquist limit using the Discrete Cosine Transform (DCT) to implicitly impose the homogeneous Neumann boundary conditions:
\begin{equation}
    \Delta u \approx \mathcal{C}^{-1} \left[ -\lambda_m \mathcal{C}[u] \right],
\end{equation}
where $\mathcal{C}$ denotes the DCT, and $\lambda_m = (m\pi/L)^2$ are the exact eigenvalues of the Laplacian operator on the discrete grid. The temporal integration is performed analytically by independently evolving each spectral mode using the exact state-transition matrix derived in Section~\ref{sec:modal_expansion}. Finally, the evolution problem of the resulting uncoupled modal equations is solved using an exact time integrator."""

tex = tex.replace(target, repl)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

