import sys

with open("Thesis.tex", "r", encoding="utf-8") as f:
    tex = f.read()

with open(r"C:\Users\gerar\.gemini\antigravity\brain\ca29bebc-1e1c-49b1-9e15-ea6a19fe6175\ard_3d_experiment.tex", "r", encoding="utf-8") as f:
    ard_3d = f.read()

target = r"\section{Differentiable Physics Neural Networks for Inverse Problems}"

if target not in tex:
    print("Target section not found!")
    sys.exit(1)

if "3D Numerical Validation" in tex:
    print("3D ARD section already exists!")
    sys.exit(0)

# Insert before the DPNN section
tex = tex.replace(target, ard_3d + "\n\n" + target)

with open("Thesis.tex", "w", encoding="utf-8") as f:
    f.write(tex)

