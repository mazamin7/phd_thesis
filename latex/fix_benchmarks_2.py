import re

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

# Fix Configuration I missing introductory sentence (if not already applied correctly)
if "\n\n\\[\nu(x,t)\n=\n\\sum_{m\\in\\{1,3\\}}\nU_m(t)\\phi_m(x),\n\\]" in tex:
    tex = tex.replace("\n\n\\[\nu(x,t)\n=\n\\sum_{m\\in\\{1,3\\}}\nU_m(t)\\phi_m(x),\n\\]", "\nThe exact solution is given by:\n\\[\nu(x,t)\n=\n\\sum_{m\\in\\{1,3\\}}\nU_m(t)\\phi_m(x).\n\\]")

# Extract and move results from Configuration V
res5 = r"This experiment evaluates the propagation accuracy of the two schemes under repeated Neumann reflections. As in Experiment~1, the FDTD method exhibits the expected exponential convergence of the discretization error. Although the Gaussian pulse is not exactly represented by a finite cosine expansion, its cosine coefficients decay exponentially. Consequently, the PSTD method rapidly reaches machine precision."
if res5 in tex:
    tex = tex.replace(res5, "")
    idx = tex.find(r"\subsection{Smooth Propagating Waves (Undamped)}")
    if idx != -1:
        # insert after the \subsection{...}
        match = re.search(r"\\subsection\{Smooth Propagating Waves \(Undamped\)\}", tex)
        if match:
            tex = tex[:match.end()] + "\n" + res5 + tex[match.end():]

# Extract and move results from Configuration VI
res6 = r"This benchmark evaluates the behavior of the two methods for propagating solutions with limited spatial regularity and sharp gradients. The FDTD method again exhibits the expected exponential convergence of the discretization error. Since the solution is non-smooth, the cosine coefficients decay only algebraically, leading the PSTD method to exhibit algebraic convergence."
if res6 in tex:
    tex = tex.replace(res6, "")
    idx = tex.find(r"\subsection{Non-smooth Propagating Waves (Undamped)}")
    if idx != -1:
        match = re.search(r"\\subsection\{Non-smooth Propagating Waves \(Undamped\)\}", tex)
        if match:
            tex = tex[:match.end()] + "\n" + res6 + tex[match.end():]
            
# Also fix comma at end of Configuration V
# The last equation in Config V is:
# \[
# x=2mL\pm\mu,
# \qquad
# m\in\mathbb Z.
# \]
# It is correct.

# Also fix the comma at the end of Configuration VI:
# Actually wait, let me check the very end of Config VI in Thesis.tex
# I need to find if there are any lingering commas.

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

