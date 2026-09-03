import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

# Add boundary conditions and zero initial velocity to the common section
target_common = r"yielding a Signal-to-Noise Ratio (SNR) of 30\,dB."
repl_common = target_common + "\nFor all experiments, we impose homogeneous Neumann boundary conditions at both ends of the domain, $\partial_x u(0,t) = \partial_x u(L,t) = 0$, and assume a zero initial velocity field, $\partial_t u(x,0) = 0$, with no external forcing, $f(x,t) = 0$."
if target_common in tex and "For all experiments, we impose homogeneous Neumann" not in tex:
    tex = tex.replace(target_common, repl_common)

# Modify Experiment 1
target_exp1 = r"\subsubsection{Experiment 1: Standing Waves}"
repl_exp1 = r"""\subsubsection{Experiment 1: Standing Waves}
In the first configuration, the initial condition is defined as a superposition of the first five cosine modes,
\[
u(x,0) = \sum_{m=1}^{5} \cos\left(\frac{m \pi x}{L}\right).
\]
"""
if target_exp1 in tex and "superposition of the first five cosine modes" not in tex:
    tex = tex.replace(target_exp1, repl_exp1)

# Modify Experiment 2
target_exp2 = r"\subsubsection{Experiment 2: Smooth Pulses}"
repl_exp2 = r"""\subsubsection{Experiment 2: Smooth Pulses}
In the second configuration, the initial condition is a smooth Gaussian pulse centered at $\mu = L/2$ with width parameter $\sigma = L/20$. To exactly satisfy the homogeneous Neumann boundary conditions without truncation artifacts, it is constructed as an even $2L$-periodic extension,
\[
u(x,0) = \sum_{m=-5}^{5} \frac{1}{2} \exp\left(-\frac{1}{2}\left(\frac{x - (2mL + \mu)}{\sigma}\right)^2\right) + \frac{1}{2} \exp\left(-\frac{1}{2}\left(\frac{x - (2mL - \mu)}{\sigma}\right)^2\right).
\]
"""
if target_exp2 in tex and "smooth Gaussian pulse centered at" not in tex:
    tex = tex.replace(target_exp2, repl_exp2)

# Modify Experiment 3
target_exp3 = r"\subsubsection{Experiment 3: Triangle Pulses (Non-smooth)}"
repl_exp3 = r"""\subsubsection{Experiment 3: Triangle Pulses (Non-smooth)}
The final experiment considers a non-smooth triangular pulse, also centered at $\mu = L/2$ with support determined by $\sigma = L/20$. Similar to the Gaussian pulse, it is constructed using an even $2L$-periodic extension,
\[
u(x,0) = \sum_{m=-5}^{5} \frac{1}{2} \max\left(1 - \frac{|x - (2mL + \mu)|}{\sigma}, 0\right) + \frac{1}{2} \max\left(1 - \frac{|x - (2mL - \mu)|}{\sigma}, 0\right).
\]
"""
if target_exp3 in tex and "non-smooth triangular pulse" not in tex:
    tex = tex.replace(target_exp3, repl_exp3)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

