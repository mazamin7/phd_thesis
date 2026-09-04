import sys

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

target1 = r"""yielding the following physical subdomains:
\begin{itemize}
    \item \textbf{Subdomain 1:} Offset $(0,0,0)$ with $30 \times 10 \times 10$ points. \\ (Physical size: $6.0 \times 2.0 \times 2.0$\,m).
    \item \textbf{Subdomain 2:} Offset $(2, 10, 0)\Delta h$ with $26 \times 10 \times 10$ points. \\ (Physical offset: $0.4 \times 2.0 \times 0.0$\,m, size: $5.2 \times 2.0 \times 2.0$\,m).
    \item \textbf{Subdomain 3:} Offset $(5, 20, 0)\Delta h$ with $20 \times 10 \times 10$ points. \\ (Physical offset: $1.0 \times 4.0 \times 0.0$\,m, size: $4.0 \times 2.0 \times 2.0$\,m).
\end{itemize}
This geometric configuration creates a sequence of narrowing connected corridors, mimicking a complex architectural layout."""

replacement1 = r"""yielding the sequence of connected subdomains depicted in Figure~\ref{fig:ard_3d_geometry}. 
\begin{figure}[H]
    \centering
    \begin{tikzpicture}[scale=1.5,
        sd/.style={draw, thick, fill=blue!10},
        dim/.style={latex-latex, thin}
        ]
        
        % Axes
        \draw[->] (-0.5,0) -- (6.5,0) node[right] {$x$ (m)};
        \draw[->] (0,-0.5) -- (0,6.5) node[above] {$y$ (m)};
        
        % Subdomain 1
        \draw[sd] (0,0) rectangle (6,2);
        \node at (3,1) {Subdomain 1};
        
        % Subdomain 2
        \draw[sd] (0.4,2) rectangle (5.6,4);
        \node at (3,3) {Subdomain 2};
        
        % Subdomain 3
        \draw[sd] (1.0,4) rectangle (5.0,6);
        \node at (3,5) {Subdomain 3};
        
        % Grid or ticks
        \foreach \x in {1,2,3,4,5,6} \draw (\x,0.05) -- (\x,-0.05) node[below] {\footnotesize \x};
        \foreach \y in {2,4,6} \draw (0.05,\y) -- (-0.05,\y) node[left] {\footnotesize \y};
        \draw (0.4, 2.05) -- (0.4, 1.95) node[below] {\footnotesize 0.4};
        \draw (5.6, 2.05) -- (5.6, 1.95) node[below] {\footnotesize 5.6};
        
        % Indicate Z thickness
        \node[anchor=west] at (4.5, 5.5) {$L_z = 2.0$\,m};
        
        % Source location
        \fill[red] (3,5) circle (1.5pt) node[right, text=red, font=\footnotesize] {Source};
        
    \end{tikzpicture}
    \caption{Top-down ($x$-$y$ plane) projection of the 3D computational domain. The geometry consists of three interconnected subdomains forming a sequence of narrowing corridors, with a uniform thickness $L_z = 2.0$\,m in the $z$-direction.}
    \label{fig:ard_3d_geometry}
\end{figure}
This geometric configuration mimics a complex architectural layout."""

tex = tex.replace(target1, replacement1)

target2 = r"where $C_c = c \frac{\Delta t}{\Delta h}$ is the Courant number"
replacement2 = r"where $\mathrm{CFL} = c \frac{\Delta t}{\Delta h}$ is the Courant--Friedrichs--Lewy (CFL) number"

target2_formula = r"S(t) = A \exp\left(-\pi^2 \left( \frac{2 C_c t}{6} - 2.0 \right)^2 \right)"
replacement2_formula = r"S(t) = A \exp\left(-\pi^2 \left( \frac{2 \, \mathrm{CFL} \, t}{6} - 2.0 \right)^2 \right)"

tex = tex.replace(target2, replacement2)
tex = tex.replace(target2_formula, replacement2_formula)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

