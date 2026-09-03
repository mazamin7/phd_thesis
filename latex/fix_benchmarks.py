import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

# Fix Configuration I
target1 = r"""u(x,0)
=
\cos\!\left(\frac{\pi x}{L}\right)
+
0.3\cos\!\left(\frac{3\pi x}{L}\right),
\]


\[
u(x,t)
=
\sum_{m\in\{1,3\}}
U_m(t)\phi_m(x),
\]"""

repl1 = r"""u(x,0)
=
\cos\!\left(\frac{\pi x}{L}\right)
+
0.3\cos\!\left(\frac{3\pi x}{L}\right).
\]
The exact solution is given by:
\[
u(x,t)
=
\sum_{m\in\{1,3\}}
U_m(t)\phi_m(x).
\]"""
tex = tex.replace(target1, repl1)


# Fix Configuration II
target2 = r"""u(x,0)
=
\sum_{j=0}^{\infty}
\frac{16}{\pi^2} \frac{\cos\left((2j+1)\pi/4\right)}{(2j+1)^2}
\cos\left((2j+1)\pi x/L\right).
\]


\[
u(x,t)
=
\sum_{j=0}^{\infty}
\frac{16}{\pi^2} \frac{\cos\left((2j+1)\pi/4\right)}{(2j+1)^2}
U_{2j+1}(t)\cos\left((2j+1)\pi x/L\right),
\]"""

repl2 = r"""u(x,0)
=
\sum_{j=0}^{\infty}
\frac{16}{\pi^2} \frac{\cos\left((2j+1)\pi/4\right)}{(2j+1)^2}
\cos\left((2j+1)\pi x/L\right).
\]
The exact solution is given by:
\[
u(x,t)
=
\sum_{j=0}^{\infty}
\frac{16}{\pi^2} \frac{\cos\left((2j+1)\pi/4\right)}{(2j+1)^2}
U_{2j+1}(t)\cos\left((2j+1)\pi x/L\right).
\]"""
tex = tex.replace(target2, repl2)


# Extract and move results from Configuration V
res5 = r"This experiment evaluates the propagation accuracy of the two schemes under repeated Neumann reflections. As in Experiment~1, the FDTD method exhibits the expected exponential convergence of the discretization error. Although the Gaussian pulse is not exactly represented by a finite cosine expansion, its cosine coefficients decay exponentially. Consequently, the PSTD method rapidly reaches machine precision."
if res5 in tex:
    tex = tex.replace(res5, "")
    # Add to numerical results
    idx = tex.find(r"\subsection{Configuration V: Smooth propagating pulse (undamped)}")
    if idx != -1:
        tex = tex[:idx+len(r"\subsection{Configuration V: Smooth propagating pulse (undamped)}")] + "\n" + res5 + tex[idx+len(r"\subsection{Configuration V: Smooth propagating pulse (undamped)}")]

# Extract and move results from Configuration VI
res6 = r"This benchmark evaluates the behavior of the two methods for propagating solutions with limited spatial regularity and sharp gradients. The FDTD method again exhibits the expected exponential convergence of the discretization error. Since the solution is non-smooth, the cosine coefficients decay only algebraically, leading the PSTD method to exhibit algebraic convergence."
if res6 in tex:
    tex = tex.replace(res6, "")
    idx = tex.find(r"\subsection{Configuration VI: Non-smooth propagating pulse (undamped)}")
    if idx != -1:
        tex = tex[:idx+len(r"\subsection{Configuration VI: Non-smooth propagating pulse (undamped)}")] + "\n" + res6 + tex[idx+len(r"\subsection{Configuration VI: Non-smooth propagating pulse (undamped)}")]

with open("Thesis_fixed.tex", "w", encoding="utf-8") as f:
    f.write(tex)

