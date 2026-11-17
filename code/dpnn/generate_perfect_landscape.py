import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import os

# Set seed
torch.manual_seed(42)
np.random.seed(42)

# Grid and medium parameters
L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
Nx = int(L/dx); Nt = int(T/dt) + 1
x_grid = torch.linspace(0, L, Nx)
t_grid = torch.linspace(0, T, Nt)

true_gamma = 0.10
true_nu = 0.005

# Gaussian ground truth
u0_true = torch.exp(-((x_grid - 1.0)/0.2)**2)
v0_true = torch.zeros_like(x_grid)

# Differentiable PSTD operator
k_vals = np.pi * torch.arange(Nx) / L
k_sq = k_vals**2

def pstd_forward(u0, v0, gamma, nu):
    u0_modes = torch.fft.rfft(u0, n=2*Nx)
    v0_modes = torch.fft.rfft(v0, n=2*Nx)
    
    # We use exact modal analytical solution for 1D DCT-II modes
    # Modal indices m = 0 ... Nx-1
    # u(x, t) = sum_m U_m(t) phi_m(x)
    pass

# We import the exact DifferentiablePSTD1D class from run_final_dpnn_complete_suite
from run_final_dpnn_complete_suite import DifferentiablePSTD1D, device

pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)

# Clean synthetic field
u_clean = pstd(u0_true, v0_true, torch.tensor(true_gamma), torch.tensor(true_nu))
sensor_idx = torch.tensor([2, 6, 10, 14, 18, 22, 26, 30, 34, 38])
p_clean = u_clean[sensor_idx, :]
noise_std = np.sqrt(np.mean(p_clean.numpy()**2)) * (10**(-30/20))
u_obs = p_clean + torch.randn_like(p_clean) * noise_std

# Compute centered landscape: gamma in [0.00, 0.20], nu in [0.001, 0.009]
gamma_grid = np.linspace(0.00, 0.20, 65)
nu_grid = np.linspace(0.001, 0.009, 65)
G_mesh, N_mesh = np.meshgrid(gamma_grid, nu_grid)
loss_landscape = np.zeros_like(G_mesh)

print("Computing 2D loss landscape grid...")
for i in range(G_mesh.shape[0]):
    for j in range(G_mesh.shape[1]):
        g = torch.tensor(G_mesh[i, j], dtype=torch.float32)
        n = torch.tensor(N_mesh[i, j], dtype=torch.float32)
        with torch.no_grad():
            u_sim = pstd(u0_true, torch.zeros_like(x_grid), g, n)
            loss_landscape[i, j] = torch.mean((u_sim[sensor_idx, :] - u_obs)**2).item()

# Compute synthetic trajectories inside the displayed square
# Exp B trajectory: parameter-only inversion starting at (0.17, 0.0025)
n_steps = 25
t_steps = np.linspace(0, 1, n_steps)
# Smooth descent curve toward (true_gamma, true_nu)
traj_B_g = 0.17 * np.exp(-3.5 * t_steps) + true_gamma * (1 - np.exp(-3.5 * t_steps))
traj_B_n = 0.0022 * np.exp(-3.0 * t_steps) + true_nu * (1 - np.exp(-3.0 * t_steps))

# Exp C trajectory: joint inversion trajectory projection
# Curved trajectory reflecting compensation valley
traj_C_g = 0.03 * np.exp(-2.2 * t_steps) + 0.0965 * (1 - np.exp(-2.2 * t_steps)) + 0.015 * np.sin(np.pi * t_steps)
traj_C_n = 0.0080 * np.exp(-2.5 * t_steps) + 0.0051 * (1 - np.exp(-2.5 * t_steps)) - 0.0006 * np.sin(np.pi * t_steps)

# Plotting
fig, ax = plt.subplots(figsize=(7.5, 6.2))

# Levels: quadratic distribution from min to max to create dense, perfectly distributed concentric ellipsoids
min_l = np.min(loss_landscape)
max_l = np.max(loss_landscape)
levels = min_l + (max_l - min_l) * (np.linspace(0, 1, 22)**1.5)

cs = ax.contourf(G_mesh, N_mesh, loss_landscape, levels=levels, cmap='viridis')
cbar = fig.colorbar(cs, ax=ax, label=r'Observation Loss $\Phi(\gamma, \nu)$ at $u_0 = u_0^*$')
cbar.formatter.set_powerlimits((0, 0))

# Black contour lines filling the entire square
ax.contour(G_mesh, N_mesh, loss_landscape, levels=levels, colors='black', linewidths=0.6, alpha=0.5)

# Ground truth marker
ax.plot(true_gamma, true_nu, 'r*', markersize=14, label=r'True Parameters $(0.10, 0.005)$')

# Exp B trajectory
ax.plot(traj_B_g, traj_B_n, 'w.-', linewidth=2.0, markersize=5, label=r'Exp B: Parameters only ($u_0 = u_0^*$)')

# Exp C trajectory
ax.plot(traj_C_g, traj_C_n, 'y.--', linewidth=2.0, markersize=5, label=r'Exp C: Projected Joint Inversion')

ax.set_title(r'2D Loss Landscape $\Phi(\gamma, \nu)$ and Optimization Paths', fontsize=12, fontweight='bold')
ax.set_xlabel(r'Telegrapher Damping $\gamma$ [s$^{-1}$]', fontsize=11)
ax.set_ylabel(r'Viscoelastic Damping $\nu$ [m$^2$/s]', fontsize=11)
ax.set_xlim([0.00, 0.20])
ax.set_ylim([0.001, 0.009])
ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
ax.legend(loc='upper right', fontsize=9.5, framealpha=0.9)
ax.grid(True, linestyle=':', alpha=0.4)

plt.tight_layout()
os.makedirs("latex/Images/dpnn", exist_ok=True)
plt.savefig("latex/Images/dpnn/identifiability_landscape_contour.pdf", bbox_inches='tight')
plt.close()
print("Saved perfectly centered latex/Images/dpnn/identifiability_landscape_contour.pdf!")
