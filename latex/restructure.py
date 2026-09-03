import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

start_idx = tex.find(r"\section{Numerical Experiments}")
end_idx = tex.find(r"\section{Optimized Schwarz Waveform Relaxation}")
num_exp_section = tex[start_idx:end_idx]

parts = re.split(r"(\\subsection\{[^\}]+\})", num_exp_section)

intro = parts[0]

benchmark_configs = ""
numerical_results = intro.replace(r"\section{Numerical Experiments}", r"\section{Numerical Results}" + "\n\\label{sec:numerical_results}\n" + "In this section, we present the numerical results for the FDTD and PSTD solvers on the benchmark configurations defined in Section~\\ref{sec:benchmark_problems}.\n")
numerical_results = numerical_results.replace("To rigorously quantify the accuracy", "As defined previously, we rigorously quantify the accuracy")

benchmark_configs += "\\section{Benchmark Problems}\n\\label{sec:benchmark_problems}\n"
benchmark_configs += "In this section, we define a fixed set of six benchmark problems that are reused throughout the thesis to evaluate different numerical methods.\n\n"
benchmark_configs += "\\subsection{Common Experimental Setup}\n"
benchmark_configs += "Unless otherwise specified, the computational domain is $\\Omega=(0,L)$ and homogeneous Neumann boundary conditions are imposed, $\\partial_xu(0,t)=\\partial_xu(L,t)=0$. The initial velocity and forcing are identically zero, $\\partial_tu(x,0)=0$, $f(x,t)=0$. Simulations are run up to a final time $T=2$.\n\n"
benchmark_configs += "To rigorously quantify the accuracy of the proposed methods, we evaluate the discrete $L^2$ error at the computational nodes at the final simulation time. Let $u(x,t)$ denote the exact analytical solution and let $u_h(x,t)$ denote the discrete numerical approximation. The discretization error is defined as\n"
benchmark_configs += "\\[\nE\n=\n\\left(\n\\Delta x\n\\sum_{i=1}^N\n|u_h(x_i,T)-u(x_i,T)|^2\n\\right)^{1/2}.\n\\]\n\n"

benchmark_configs += "\\subsection{Reference Solutions}\n"
benchmark_configs += "The exact reference solution for each configuration is obtained by evolving the modal amplitudes analytically (see Section~\\ref{sec:modal_expansion}). Since the initial conditions are typically linear combinations of cosine eigenfunctions, the exact solution is obtained by evolving each mode independently. For non-smooth profiles, the cosine series is truncated after a sufficiently large number of modes, making the truncation error negligible.\n\n"

benchmark_configs += "\\subsection{Benchmark Configurations}\n"

exp_titles = [
    "Configuration I: Smooth standing wave (undamped)",
    "Configuration II: Non-smooth standing wave (undamped)",
    "Configuration III: Smooth standing wave (damped)",
    "Configuration IV: Non-smooth standing wave (damped)",
    "Configuration V: Smooth propagating pulse (undamped)",
    "Configuration VI: Non-smooth propagating pulse (undamped)"
]

markers = [
    "This benchmark verifies",
    "This experiment investigates",
    "This configuration evaluates",
    "Similar to the undamped case,",
    "This experiment models",
    "The final experiment evaluates"
]

for i in range(1, len(parts), 2):
    title = parts[i]
    content = parts[i+1]
    
    # Find the split point
    split_point = -1
    for m in markers:
        idx = content.find(m)
        if idx != -1:
            if split_point == -1 or idx < split_point:
                split_point = idx
                
    if split_point == -1:
        # If no marker is found, just dump everything before \begin{figure} into definition
        fig_idx = content.find(r"\begin{figure}")
        if fig_idx != -1:
            split_point = fig_idx
        else:
            split_point = len(content)
            
    definition = content[:split_point].strip()
    results = content[split_point:].strip()
    
    # Clean up some redundant definitions
    # Replace "The computational domain is \Omega=(0,L) and homogeneous Neumann boundary conditions are imposed," with ""
    # etc. But keeping it as is is also fine and mathematically rigorous. Let's just remove the first redundant lines if they match
    lines = definition.split('\n')
    filtered_lines = []
    skip_next = False
    for j, line in enumerate(lines):
        if "The computational domain is" in line or "The computational domain," in line:
            skip_next = True
            continue
        if skip_next and line.strip() == "\\]":
            skip_next = False
            continue
        if skip_next:
            continue
            
        if "The damping coefficients are set to" in line or "The initial velocity and forcing" in line:
            skip_next = True
            continue
            
        filtered_lines.append(line)
        
    cleaned_def = "\n".join(filtered_lines)
    
    idx = (i-1)//2
    benchmark_configs += f"\\subsubsection{{{exp_titles[idx]}}}\n"
    # we just use the raw definition for safety to not lose equations
    # but we can filter some known exact matches
    # actually, raw definition is safer! Let's just use it
    
    # Let's remove the common redundant parts manually:
    definition = definition.replace("The computational domain is $\\Omega=(0,L)$ and homogeneous Neumann boundary conditions are imposed,\n\\[\n\\partial_xu(0,t)=\\partial_xu(L,t)=0.\n\\]", "")
    definition = definition.replace("The computational domain, boundary conditions and model parameters are the same as in the previous experiment,\n\\[\n\\Omega=(0,L),\\qquad\n\\partial_xu(0,t)=\\partial_xu(L,t)=0,\\qquad\n\\gamma=\\nu=0.\n\\]", "")
    definition = definition.replace("while the initial velocity and forcing are identically zero,\n\\[\n\\partial_tu(x,0)=0,\n\\qquad\nf(x,t)=0.\n\\]", "")
    definition = definition.replace("The initial velocity and forcing vanish,\n\\[\n\\partial_tu(x,0)=0,\n\\qquad\nf(x,t)=0.\n\\]", "")
    definition = definition.replace("The reference solution is obtained by evolving the modal amplitudes,", "")
    definition = definition.replace("where the modal amplitudes $U_{2j+1}(t)$ evolve according to the exact modal solution. In practice, the series is truncated after a sufficiently large number of modes, making the truncation error negligible.", "")
    definition = definition.replace("Since the initial condition is a linear combination of two cosine eigenfunctions, the exact solution is obtained by evolving each mode independently,", "")
    definition = definition.replace("where $U_m(t)$ denotes the time evolution of mode $m$ as shown in Section~\\ref{sec:modal_expansion}.", "")
    definition = definition.strip()
    
    benchmark_configs += definition + "\n\n"
    
    numerical_results += title + "\n"
    numerical_results += results + "\n\n"

# 2. Insert benchmark_configs before \section{Finite-Difference Time-Domain Method}
fdtd_idx = tex.find(r"\section{Finite-Difference Time-Domain Method}")
tex = tex[:fdtd_idx] + benchmark_configs + "\n\n" + tex[fdtd_idx:]

# 3. Replace the original Numerical Experiments section
start_idx = tex.find(r"\section{Numerical Experiments}")
end_idx = tex.find(r"\section{Optimized Schwarz Waveform Relaxation}")
tex = tex[:start_idx] + numerical_results + tex[end_idx:]

# 4. Fix OSWR section
oswr_str = r"""\subsection{Numerical Experiments}

We employ the monolithic FDTD scheme introduced in the previous chapter. To validate the theoretical convergence rates of the Optimized Schwarz Waveform Relaxation (OSWR) method, we evaluate the $L^\infty$ domain decomposition error (defined as the maximum absolute difference between the SWR solution and the monolithic FDTD reference solution) over successive Schwarz iterations across six distinct experiments:

\begin{enumerate}
    \item Smooth standing wave (undamped)
    \item Non-smooth standing wave (undamped)
    \item Smooth standing wave (damped)
    \item Non-smooth standing wave (damped)
    \item Smooth propagating pulse (undamped)
    \item Non-smooth propagating pulse (undamped)
\end{enumerate}

The domain $\Omega=(0,L)$ is split into two overlapping subdomains with an overlap size $\delta=L/10$. The physical parameters vary across the configurations, including highly damped and undamped cases. The domain is discretized using a finite-difference spatial step $\Delta x=0.05$ and time step $\Delta t=0.05$. In all configurations, the simulations are run for a total time duration $T=2$. It is important to note that for spacetime domain decomposition methods like SWR, the total number of iterations required to reach convergence inherently depends on this simulation duration $T$, as boundary information must propagate back and forth across the entire spacetime interface."""

oswr_repl = r"""\subsection{Numerical Results}

We employ the monolithic FDTD scheme introduced in the previous chapter. To validate the theoretical convergence rates of the Optimized Schwarz Waveform Relaxation (OSWR) method, we evaluate the $L^\infty$ domain decomposition error (defined as the maximum absolute difference between the SWR solution and the monolithic FDTD reference solution) over successive Schwarz iterations across the six benchmark configurations defined in Section~\ref{sec:benchmark_problems}.

The domain $\Omega=(0,L)$ is split into two overlapping subdomains with an overlap size $\delta=L/10$. The physical parameters, grid resolution ($\Delta x=0.05$, $\Delta t=0.05$), and total simulation duration ($T=2$) are identical to those defined in the baseline configurations. It is important to note that for spacetime domain decomposition methods like SWR, the total number of iterations required to reach convergence inherently depends on this simulation duration $T$, as boundary information must propagate back and forth across the entire spacetime interface."""

tex = tex.replace(oswr_str, oswr_repl)

# 5. Fix ARD section
ard_str = r"""as a function of the interface spatial order and grid resolution. The domain $\Omega=(0,L)$ is decomposed into two non-overlapping subdomains, and the ARD interface operator is applied with varying stencil sizes (Order 2, 4, 6, and 8). The experiments encompass the same six configurations with differing source smoothness and damping profiles utilized in the previous section:

\begin{enumerate}
    \item Smooth standing wave (undamped)
    \item Non-smooth standing wave (undamped)
    \item Smooth standing wave (damped)
    \item Non-smooth standing wave (damped)
    \item Smooth propagating pulse (undamped)
    \item Non-smooth propagating pulse (undamped)
\end{enumerate}

For the propagating pulses (Experiments 5 and 6), they are placed asymmetrically at $x=L/4$ to prevent artificial symmetry at the boundary."""

ard_repl = r"""as a function of the interface spatial order and grid resolution. The domain $\Omega=(0,L)$ is decomposed into two non-overlapping subdomains, and the ARD interface operator is applied with varying stencil sizes (Order 2, 4, 6, and 8). We evaluate the method across the six benchmark configurations defined in Section~\ref{sec:benchmark_problems}. For the propagating pulses (Configurations V and VI), they are placed asymmetrically at $x=L/4$ to prevent artificial symmetry at the boundary."""

tex = tex.replace(ard_str, ard_repl)

with open("Thesis_restructured.tex", "w", encoding="utf-8") as f:
    f.write(tex)

