import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

res5 = r"This experiment evaluates the propagation accuracy of the two schemes under repeated Neumann reflections. As in Experiment~1, the FDTD method exhibits the expected exponential convergence of the discretization error. Although the Gaussian pulse is not exactly represented by a finite cosine expansion, its cosine coefficients decay exponentially. Consequently, the PSTD method rapidly reaches machine precision."
target5 = r"\subsection{Smooth Propagating Waves (Undamped)}"
if target5 in tex and res5 not in tex:
    tex = tex.replace(target5, target5 + "\n" + res5)

res6 = r"This benchmark evaluates the behavior of the two methods for propagating solutions with limited spatial regularity and sharp gradients. The FDTD method again exhibits the expected exponential convergence of the discretization error. Since the solution is non-smooth, the cosine coefficients decay only algebraically, leading the PSTD method to exhibit algebraic convergence."
target6 = r"\subsection{Non-smooth Propagating Waves (Undamped)}"
if target6 in tex and res6 not in tex:
    tex = tex.replace(target6, target6 + "\n" + res6)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

