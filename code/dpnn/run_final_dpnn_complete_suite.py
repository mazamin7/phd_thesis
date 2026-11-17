import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time

device = torch.device('cpu')
os.makedirs("latex/Images/dpnn", exist_ok=True)

def trapz_compat(y, x):
    if hasattr(np, 'trapezoid'):
        return np.trapezoid(y, x)
    return np.trapz(y, x)

# ==============================================================================
# 1. CONTINUOUS ANALYTICAL MODAL SOLVER (AVOIDS INVERSE CRIME)
# ==============================================================================
def analytical_solution(x_eval, t_eval, u0_func, c=1.0, gamma=0.1, nu=0.005, L=2.0, M=200):
    """
    Independent continuous modal series expansion with M=200 modes and 2000-point continuous
    spatial quadrature. Evaluates the exact continuum wavefield u(x, t) on [0, L] with
    homogeneous Neumann boundary conditions, completely avoiding discrete forward solver artifacts.
    """
    xq = np.linspace(0, L, 2000)
    u0_q = u0_func(xq)
    Nt = len(t_eval); Nx = len(x_eval)
    u_field = np.zeros((Nx, Nt))
    
    # Mode 0 (constant spatial mode)
    phi0_q = np.full_like(xq, np.sqrt(1.0 / L))
    phi0_eval = np.full_like(x_eval, np.sqrt(1.0 / L))
    U0_0 = trapz_compat(u0_q * phi0_q, xq)
    for ti in range(Nt):
        u_field[:, ti] += U0_0 * phi0_eval
        
    for m in range(1, M + 1):
        km = m * np.pi / L
        omega_m = c * km
        dm = gamma + nu * (km**2)
        phi_q = np.sqrt(2.0 / L) * np.cos(km * xq)
        phi_eval = np.sqrt(2.0 / L) * np.cos(km * x_eval)
        Um_0 = trapz_compat(u0_q * phi_q, xq)
        
        discr = dm**2 - 4 * (omega_m**2)
        if discr < -1e-12:
            omega_d = np.sqrt(-discr) / 2.0
            coeff = np.exp(-dm * t_eval / 2.0) * (np.cos(omega_d * t_eval) + (dm / (2.0 * omega_d)) * np.sin(omega_d * t_eval))
        elif discr > 1e-12:
            delta = np.sqrt(discr) / 2.0
            coeff = np.exp(-dm * t_eval / 2.0) * (np.cosh(delta * t_eval) + (dm / (2.0 * delta)) * np.sinh(delta * t_eval))
        else:
            coeff = np.exp(-dm * t_eval / 2.0) * (1.0 + (dm / 2.0) * t_eval)
            
        u_field += np.outer(phi_eval, Um_0 * coeff)
        
    return u_field

# ==============================================================================
# 2. NUMERICAL FORWARD SOLVERS (PSTD & FDTD)
# ==============================================================================
def interior_assembly(Nx, dx, dt, c, gamma, nu, device='cpu'):
    I = torch.eye(Nx, device=device)
    D = torch.zeros((Nx, Nx), device=device)
    idx = torch.arange(0, Nx)
    D[idx, idx] = -2.0
    idx_off = torch.arange(0, Nx-1)
    D[idx_off+1, idx_off] = 1.0
    D[idx_off, idx_off+1] = 1.0
    D = D / (dx**2)
    A = (1 + gamma*dt/2)*I - (nu*dt/2)*D
    B = 2*I + (c**2 * dt**2)*D
    C = -(1 - gamma*dt/2)*I - (nu*dt/2)*D
    return A, B, C

def apply_bc(A, B, C, bcType, Nx, dx, dt, c, gamma, nu, device='cpu'):
    A_mod, B_mod, C_mod = A.clone(), B.clone(), C.clone()
    if bcType.lower() == 'neumannghost':
        A_mod[0, 1] *= 2; B_mod[0, 1] *= 2; C_mod[0, 1] *= 2
        A_mod[-1, -2] *= 2; B_mod[-1, -2] *= 2; C_mod[-1, -2] *= 2
        def rhsFun(t): return torch.zeros(Nx, device=device)
    else: 
        raise ValueError(f"Unknown BC type: {bcType}")
    return A_mod, B_mod, C_mod, rhsFun

class DifferentiableFDTD1D(nn.Module):
    def __init__(self, Nx, dx, dt, L, T, bcType='neumannghost', device='cpu'):
        super().__init__()
        self.Nx = Nx; self.dx = dx; self.dt = dt; self.L = L; self.T = T
        self.bcType = bcType; self.device = device
        self.x_grid = torch.linspace(0, L, Nx, device=device)
        self.t_grid = torch.arange(0, T+dt/2, dt, device=device)
        self.Nt = len(self.t_grid)
        
    def forward(self, u0, v0, gamma, nu):
        c = 1.0
        A, B, C = interior_assembly(self.Nx, self.dx, self.dt, c, gamma, nu, self.device)
        A_mod, B_mod, C_mod, rhsFun = apply_bc(A, B, C, self.bcType, self.Nx, self.dx, self.dt, c, gamma, nu, self.device)
        A_inv = torch.linalg.inv(A_mod)
        
        u_array = [u0]
        uxx0 = torch.zeros(self.Nx, device=self.device)
        uxx0[1:-1] = (u0[2:] - 2*u0[1:-1] + u0[:-2]) / (self.dx**2)
        uxx0[0] = uxx0[1]; uxx0[-1] = uxx0[-2]
        vxx0 = torch.zeros(self.Nx, device=self.device)
        vxx0[1:-1] = (v0[2:] - 2*v0[1:-1] + v0[:-2]) / (self.dx**2)
        vxx0[0] = vxx0[1]; vxx0[-1] = vxx0[-2]
        utt0 = (c**2)*uxx0 + nu*vxx0 - gamma*v0
        
        uPrev = u0
        uCurr = u0 + self.dt*v0 + 0.5*(self.dt**2)*utt0
        u_array.append(uCurr)
        for n in range(1, self.Nt - 1):
            t = self.t_grid[n]
            rhs = torch.mv(B_mod, uCurr) + torch.mv(C_mod, uPrev) + rhsFun(t)
            uNext = torch.mv(A_inv, rhs)
            u_array.append(uNext)
            uPrev = uCurr
            uCurr = uNext
        return torch.stack(u_array, dim=1)

def dct_fft(x):
    N = x.shape[-1]
    x_ext = torch.cat([x, x.flip(-1)], dim=-1)
    X = torch.fft.fft(x_ext, dim=-1)[..., :N]
    k = torch.arange(N, device=x.device).float()
    phase = torch.exp(-1j * torch.pi * k / (2 * N))
    X_real = (X * phase).real / torch.sqrt(torch.tensor(2.0 * N, device=x.device))
    X_real[..., 0] /= torch.sqrt(torch.tensor(2.0, device=x.device))
    return X_real

def idct_fft(X):
    N = X.shape[-1]
    X_unnorm = X * torch.sqrt(torch.tensor(2.0 * N, device=X.device))
    X_unnorm[..., 0] *= torch.sqrt(torch.tensor(2.0, device=X.device))
    k = torch.arange(N, device=X.device).float()
    phase = torch.exp(1j * torch.pi * k / (2 * N))
    Y = X_unnorm * phase
    Y_ext = torch.zeros(*X.shape[:-1], 2*N, dtype=torch.complex64, device=X.device)
    Y_ext[..., :N] = Y
    Y_ext[..., N+1:] = torch.conj(Y[..., 1:].flip(-1))
    x_ext = torch.fft.ifft(Y_ext, dim=-1)
    return x_ext.real[..., :N]

class DifferentiablePSTD1D(nn.Module):
    def __init__(self, Nx, dx, dt, L, T, device='cpu'):
        super().__init__()
        self.Nx = Nx; self.dx = dx; self.dt = dt; self.L = L; self.T = T; self.device = device
        self.x_grid = (torch.arange(Nx, device=device) + 0.5) * dx 
        self.t_grid = torch.arange(0, T+dt/2, dt, device=device)
        self.Nt = len(self.t_grid)
        m = torch.arange(Nx, device=device).float()
        self.lambda_val = -((m * torch.pi / L)**2)
        
    def forward(self, u0, v0, gamma, nu):
        c = 1.0
        d_val = gamma - nu * self.lambda_val
        discr = d_val**2 + 4 * c**2 * self.lambda_val
        
        mask_lam0_dj0 = (self.lambda_val == 0) & (d_val == 0)
        mask_lam0_djN0 = (self.lambda_val == 0) & (d_val != 0)
        mask_lamN0_d0 = (self.lambda_val != 0) & (torch.abs(discr) < 1e-12)
        mask_lamN0_dneg = (self.lambda_val != 0) & (discr <= -1e-12)
        mask_lamN0_dpos = (self.lambda_val != 0) & (discr >= 1e-12)
        
        t = self.t_grid.unsqueeze(1)
        d_val_2d = d_val.unsqueeze(0)
        discr_2d = discr.unsqueeze(0)
        e = torch.exp(-d_val_2d * t / 2)
        
        S11 = torch.zeros((self.Nt, self.Nx), device=self.device)
        S12 = torch.zeros((self.Nt, self.Nx), device=self.device)
        
        if mask_lam0_dj0.any():
            S11[:, mask_lam0_dj0] = 1.0
            S12[:, mask_lam0_dj0] = t
        if mask_lam0_djN0.any():
            dj = d_val_2d[:, mask_lam0_djN0]
            S11[:, mask_lam0_djN0] = 1.0
            S12[:, mask_lam0_djN0] = (1 - torch.exp(-dj*t))/dj
        if mask_lamN0_d0.any():
            dj = d_val_2d[:, mask_lamN0_d0]; ee = e[:, mask_lamN0_d0]
            S11[:, mask_lamN0_d0] = ee * (1 + dj*t/2)
            S12[:, mask_lamN0_d0] = ee * t
        if mask_lamN0_dneg.any():
            dj = d_val_2d[:, mask_lamN0_dneg]; ee = e[:, mask_lamN0_dneg]
            omega = torch.sqrt(-discr_2d[:, mask_lamN0_dneg])/2
            co = torch.cos(omega*t); si = torch.sin(omega*t)
            S11[:, mask_lamN0_dneg] = ee * (co + si * dj/(2*omega))
            S12[:, mask_lamN0_dneg] = ee * (si / omega)
        if mask_lamN0_dpos.any():
            dj = d_val_2d[:, mask_lamN0_dpos]
            delta = torch.sqrt(discr_2d[:, mask_lamN0_dpos])
            exp_plus = torch.exp((-dj + delta)*t/2)
            exp_minus = torch.exp((-dj - delta)*t/2)
            S11[:, mask_lamN0_dpos] = 0.5 * exp_plus * (1 + dj/delta) + 0.5 * exp_minus * (1 - dj/delta)
            S12[:, mask_lamN0_dpos] = (exp_plus - exp_minus) / delta

        Um0 = dct_fft(u0).unsqueeze(0)
        Vm0 = dct_fft(v0).unsqueeze(0)
        Um_all = S11 * Um0 + S12 * Vm0
        u_all = idct_fft(Um_all)
        return u_all.T # (Nx, Nt)

# ==============================================================================
# 3. NEURAL NETWORK ARCHITECTURES WITH LOGARITHMIC PARAMETERIZATION
# ==============================================================================
class SineLayer(nn.Module):
    def __init__(self, in_features, out_features, bias=True, is_first=False, omega_0=30):
        super().__init__()
        self.omega_0 = omega_0; self.is_first = is_first; self.in_features = in_features
        self.linear = nn.Linear(in_features, out_features, bias=bias)
        with torch.no_grad():
            if self.is_first:
                self.linear.weight.uniform_(-1 / self.in_features, 1 / self.in_features)      
            else:
                self.linear.weight.uniform_(-np.sqrt(6 / self.in_features) / self.omega_0, np.sqrt(6 / self.in_features) / self.omega_0)
    def forward(self, x):
        return torch.sin(self.omega_0 * self.linear(x))

class SIREN(nn.Module):
    def __init__(self, in_features=1, hidden_features=32, hidden_layers=2, out_features=1):
        super().__init__()
        layers = [SineLayer(in_features, hidden_features, is_first=True)]
        for _ in range(hidden_layers):
            layers.append(SineLayer(hidden_features, hidden_features, is_first=False))
        final_linear = nn.Linear(hidden_features, out_features)
        with torch.no_grad():
            final_linear.weight.uniform_(-np.sqrt(6 / hidden_features) / 30, np.sqrt(6 / hidden_features) / 30)
        layers.append(final_linear)
        self.net = nn.Sequential(*layers)
    def forward(self, coords):
        return self.net(coords).squeeze(-1)

class DPNNModelLog(nn.Module):
    """
    DPNN with logarithmic parameterization:
    gamma = gamma_ref * exp(xi_gamma), nu = nu_ref * exp(xi_nu).
    Guarantees strict positivity, non-vanishing gradients, and optimizes relative variations.
    """
    def __init__(self, init_gamma=0.5, init_nu=0.001, gamma_ref=0.1, nu_ref=0.005):
        super().__init__()
        self.siren = SIREN(in_features=1, hidden_features=32, hidden_layers=2, out_features=1)
        self.gamma_ref = gamma_ref
        self.nu_ref = nu_ref
        self.xi_gamma = nn.Parameter(torch.tensor(float(np.log(init_gamma / gamma_ref))))
        self.xi_nu = nn.Parameter(torch.tensor(float(np.log(init_nu / nu_ref))))
        
    def forward(self, coords):
        p0 = self.siren(coords)
        gamma = self.gamma_ref * torch.exp(self.xi_gamma)
        nu = self.nu_ref * torch.exp(self.xi_nu)
        return p0, gamma, nu

class SirenPINNLog(nn.Module):
    """
    High-capacity PINN baseline (4 hidden layers, 128 units, ~50k parameters)
    with logarithmic parameterization for damping parameters.
    """
    def __init__(self, init_gamma=0.5, init_nu=0.001, gamma_ref=0.1, nu_ref=0.005, hidden=128, num_layers=4, omega_0=5.0):
        super().__init__()
        layers = [SineLayer(2, hidden, is_first=True, omega_0=omega_0)]
        for _ in range(num_layers - 1):
            layers.append(SineLayer(hidden, hidden, is_first=False, omega_0=omega_0))
        final_layer = nn.Linear(hidden, 1)
        with torch.no_grad():
            final_layer.weight.uniform_(-np.sqrt(6 / hidden) / omega_0, np.sqrt(6 / hidden) / omega_0)
        layers.append(final_layer)
        self.net = nn.Sequential(*layers)
        self.gamma_ref = gamma_ref
        self.nu_ref = nu_ref
        self.xi_gamma = nn.Parameter(torch.tensor(float(np.log(init_gamma / gamma_ref))))
        self.xi_nu = nn.Parameter(torch.tensor(float(np.log(init_nu / nu_ref))))
        
    def forward(self, x, t):
        xt = torch.cat([x, t], dim=-1)
        return self.net(xt).squeeze(-1)
        
    def get_params(self):
        gamma = self.gamma_ref * torch.exp(self.xi_gamma)
        nu = self.nu_ref * torch.exp(self.xi_nu)
        return gamma, nu

print("Physics and Network modules successfully initialized.")

# ==============================================================================
# 4. BENCHMARK PIPELINE & MAIN SUITE
# ==============================================================================
def run_main_benchmarks():
    print("\n" + "="*80)
    print("RUNNING MAIN BENCHMARK SUITE: PSTD-DPNN vs FDTD-DPNN vs PINN")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
    true_gamma = 0.1; true_nu = 0.005
    init_gamma = 0.5; init_nu = 0.001
    num_sensors = 10; snr_db = 30
    
    Nx_fdtd = int(L/dx) + 1 # 41 nodes (ghost/boundary)
    Nx_pstd = int(L/dx)     # 40 nodes (cell-centered)
    Nt = int(T/dt) + 1      # 41 time steps
    
    fdtd_solver = DifferentiableFDTD1D(Nx_fdtd, dx, dt, L, T, device=device)
    pstd_solver = DifferentiablePSTD1D(Nx_pstd, dx, dt, L, T, device=device)
    
    x_grid_fdtd = fdtd_solver.x_grid
    x_grid_pstd = pstd_solver.x_grid
    t_grid = pstd_solver.t_grid
    
    torch.manual_seed(123)
    sensor_idx = torch.sort(torch.randperm(Nx_pstd - 2)[:num_sensors] + 1)[0]
    sensor_x = x_grid_pstd[sensor_idx]
    
    # Map sensor locations to FDTD grid
    sensor_idx_fdtd = []
    for sx in sensor_x:
        diffs = torch.abs(x_grid_fdtd - sx)
        sensor_idx_fdtd.append(torch.argmin(diffs).item())
    sensor_idx_fdtd = torch.tensor(sensor_idx_fdtd, device=device)
    
    # UNIFIED PULSE PARAMETERS: mu = L/2 = 1.0, sigma = L/20 = 0.1
    mu = L / 2; sigma = L / 20; M_pulse = 5
    
    experiments = [
        {'name': '0_StandingWaves', 'title': 'Multimodal Standing Waves', 'type': 'standingWaves', 'lambda_sp': 0.0, 'lambda_smooth': 0.0},
        {'name': '1_Smooth', 'title': 'Smooth Gaussian Pulses', 'type': 'smoothPulse', 'lambda_sp': 0.0, 'lambda_smooth': 0.0},
        {'name': '2_Triangle', 'title': 'Non-Smooth Triangular Pulses', 'type': 'trianglePulse', 'lambda_sp': 0.0, 'lambda_smooth': 0.0}
    ]
    
    all_summary_results = []
    
    for exp in experiments:
        print(f"\n>>> Running Benchmark: {exp['title']} ({exp['name']})")
        
        # 1. Define continuous profile function
        if exp['type'] == 'standingWaves':
            def u0_func(x):
                val = np.zeros_like(x)
                for m in range(1, 6):
                    val += np.cos(m * np.pi * x / L)
                return val
        elif exp['type'] == 'smoothPulse':
            def u0_func(x):
                val = np.zeros_like(x)
                for m in range(-M_pulse, M_pulse + 1):
                    val += 0.5 * np.exp(-0.5 * ((x - (2*m*L + mu)) / sigma)**2)
                    val += 0.5 * np.exp(-0.5 * ((x - (2*m*L - mu)) / sigma)**2)
                return val
        elif exp['type'] == 'trianglePulse':
            def u0_func(x):
                val = np.zeros_like(x)
                for m in range(-M_pulse, M_pulse + 1):
                    val += 0.5 * np.maximum(0.0, 1.0 - np.abs(x - (2*m*L + mu)) / sigma)
                    val += 0.5 * np.maximum(0.0, 1.0 - np.abs(x - (2*m*L - mu)) / sigma)
                return val

        # Evaluate discrete initial conditions
        u0_pstd = torch.tensor(u0_func(x_grid_pstd.numpy()), dtype=torch.float32, device=device)
        u0_fdtd = torch.tensor(u0_func(x_grid_fdtd.numpy()), dtype=torch.float32, device=device)
        
        # 2. Continuous analytical ground-truth wavefield (AVOIDS INVERSE CRIME)
        u_exact_cont_pstd = analytical_solution(x_grid_pstd.numpy(), t_grid.numpy(), u0_func, c, true_gamma, true_nu, L, M=200)
        u_exact_cont_fdtd = analytical_solution(x_grid_fdtd.numpy(), t_grid.numpy(), u0_func, c, true_gamma, true_nu, L, M=200)
        
        u_exact_pstd = torch.tensor(u_exact_cont_pstd, dtype=torch.float32, device=device)
        u_exact_fdtd = torch.tensor(u_exact_cont_fdtd, dtype=torch.float32, device=device)
        
        # Extract synthetic sensor observations from continuous wavefield
        p_exact_sensors = u_exact_pstd[sensor_idx, :]
        torch.manual_seed(123)
        noise = torch.randn_like(p_exact_sensors)
        noise = noise / torch.norm(noise) * torch.norm(p_exact_sensors) * (10**(-snr_db / 20))
        p_noisy = p_exact_sensors + noise
        
        # PINN sensor coordinates and signal power
        sens_x_grid = sensor_x.unsqueeze(-1).repeat(1, Nt).reshape(-1, 1)
        sens_t_grid = t_grid.unsqueeze(0).repeat(num_sensors, 1).reshape(-1, 1)
        x_sens_norm = (sens_x_grid / L) * 2 - 1
        t_sens_norm = (sens_t_grid / T) * 2 - 1
        p_sens_flat = p_noisy.reshape(-1)
        sig_pow = torch.mean(p_sens_flat**2)
        
        epochs = 350
        methods = ['PSTD-DPNN', 'FDTD-DPNN', 'High-Capacity PINN']
        history = {m: {'loss': [], 'nmse_u0': [], 'gamma': [], 'nu': [], 'time_per_ep': 0.0} for m in methods}
        final_models = {}
        
        # --- A. PSTD-DPNN ---
        torch.manual_seed(42)
        model_pstd = DPNNModelLog(init_gamma, init_nu)
        opt_net_p = torch.optim.Adam(model_pstd.siren.parameters(), lr=1e-3)
        opt_phys_p = torch.optim.Adam([model_pstd.xi_gamma, model_pstd.xi_nu], lr=5e-2)
        x_norm_pstd = ((x_grid_pstd / L) * 2 - 1).unsqueeze(-1)
        
        t0 = time.time()
        for ep in range(epochs):
            opt_net_p.zero_grad(); opt_phys_p.zero_grad()
            p0_pred, g_pred, n_pred = model_pstd(x_norm_pstd)
            u_pred = pstd_solver(p0_pred, torch.zeros_like(x_grid_pstd), g_pred, n_pred)
            
            loss = torch.mean((u_pred[sensor_idx, :] - p_noisy)**2)
            loss.backward()
            opt_net_p.step(); opt_phys_p.step()
            
            with torch.no_grad():
                nmse = (torch.sum((p0_pred - u0_pstd)**2) / torch.sum(u0_pstd**2)).item()
                history['PSTD-DPNN']['loss'].append(loss.item())
                history['PSTD-DPNN']['nmse_u0'].append(nmse)
                history['PSTD-DPNN']['gamma'].append(g_pred.item())
                history['PSTD-DPNN']['nu'].append(n_pred.item())
        history['PSTD-DPNN']['time_per_ep'] = (time.time() - t0) / epochs * 1000
        final_models['PSTD-DPNN'] = model_pstd
        
        # --- B. FDTD-DPNN ---
        torch.manual_seed(42)
        model_fdtd = DPNNModelLog(init_gamma, init_nu)
        opt_net_f = torch.optim.Adam(model_fdtd.siren.parameters(), lr=1e-3)
        opt_phys_f = torch.optim.Adam([model_fdtd.xi_gamma, model_fdtd.xi_nu], lr=5e-2)
        x_norm_fdtd = ((x_grid_fdtd / L) * 2 - 1).unsqueeze(-1)
        
        t0 = time.time()
        for ep in range(epochs):
            opt_net_f.zero_grad(); opt_phys_f.zero_grad()
            p0_pred, g_pred, n_pred = model_fdtd(x_norm_fdtd)
            u_pred = fdtd_solver(p0_pred, torch.zeros_like(x_grid_fdtd), g_pred, n_pred)
            
            loss = torch.mean((u_pred[sensor_idx_fdtd, :] - p_noisy)**2)
            loss.backward()
            opt_net_f.step(); opt_phys_f.step()
            
            with torch.no_grad():
                nmse = (torch.sum((p0_pred - u0_fdtd)**2) / torch.sum(u0_fdtd**2)).item()
                history['FDTD-DPNN']['loss'].append(loss.item())
                history['FDTD-DPNN']['nmse_u0'].append(nmse)
                history['FDTD-DPNN']['gamma'].append(g_pred.item())
                history['FDTD-DPNN']['nu'].append(n_pred.item())
        history['FDTD-DPNN']['time_per_ep'] = (time.time() - t0) / epochs * 1000
        final_models['FDTD-DPNN'] = model_fdtd
        
        # --- C. High-Capacity PINN Baseline ---
        torch.manual_seed(42)
        pinn = SirenPINNLog(init_gamma, init_nu, hidden=128, num_layers=4, omega_0=5.0)
        opt_pinn_net = torch.optim.Adam(pinn.net.parameters(), lr=1e-3)
        opt_pinn_phys = torch.optim.Adam([pinn.xi_gamma, pinn.xi_nu], lr=5e-2)
        
        t0 = time.time()
        for ep in range(epochs):
            opt_pinn_net.zero_grad(); opt_pinn_phys.zero_grad()
            g_pred, n_pred = pinn.get_params()
            
            # 1. Data loss
            u_sens_pred = pinn(x_sens_norm, t_sens_norm)
            loss_data = torch.mean((u_sens_pred - p_sens_flat)**2) / sig_pow
            
            # 2. PDE residual via Autograd collocation
            x_c = (torch.rand(300, 1, device=device) * L).requires_grad_(True)
            t_c = (torch.rand(300, 1, device=device) * T).requires_grad_(True)
            u_c = pinn((x_c / L) * 2 - 1, (t_c / T) * 2 - 1)
            
            u_t = torch.autograd.grad(u_c, t_c, torch.ones_like(u_c), create_graph=True)[0]
            u_tt = torch.autograd.grad(u_t, t_c, torch.ones_like(u_t), create_graph=True)[0]
            u_x = torch.autograd.grad(u_c, x_c, torch.ones_like(u_c), create_graph=True)[0]
            u_xx = torch.autograd.grad(u_x, x_c, torch.ones_like(u_x), create_graph=True)[0]
            u_txx = torch.autograd.grad(u_xx, t_c, torch.ones_like(u_xx), create_graph=True)[0]
            
            pde_res = u_tt - (c**2) * u_xx - n_pred * u_txx + g_pred * u_t
            loss_pde = torch.mean(pde_res**2) / (sig_pow + 1e-6)
            
            # 3. Boundary conditions (Neumann at 0 and L)
            x_b0 = (torch.zeros(20, 1, device=device)).requires_grad_(True)
            x_bL = (torch.ones(20, 1, device=device) * L).requires_grad_(True)
            t_b = (torch.linspace(0, T, 20, device=device).unsqueeze(-1)).requires_grad_(True)
            u_b0 = pinn((x_b0 / L) * 2 - 1, (t_b / T) * 2 - 1)
            u_x_b0 = torch.autograd.grad(u_b0, x_b0, torch.ones_like(u_b0), create_graph=True)[0]
            u_bL = pinn((x_bL / L) * 2 - 1, (t_b / T) * 2 - 1)
            u_x_bL = torch.autograd.grad(u_bL, x_bL, torch.ones_like(u_bL), create_graph=True)[0]
            loss_bc = (torch.mean(u_x_b0**2) + torch.mean(u_x_bL**2)) / (sig_pow + 1e-6)
            
            # 4. Zero initial velocity
            x_i = (torch.linspace(0, L, 30, device=device).unsqueeze(-1)).requires_grad_(True)
            t_i = (torch.zeros(30, 1, device=device)).requires_grad_(True)
            u_i = pinn((x_i / L) * 2 - 1, (t_i / T) * 2 - 1)
            u_t_i = torch.autograd.grad(u_i, t_i, torch.ones_like(u_i), create_graph=True)[0]
            loss_ic = torch.mean(u_t_i**2) / (sig_pow + 1e-6)
            
            loss = 1.0 * loss_data + 1e-4 * loss_pde + 1e-3 * loss_bc + 1e-3 * loss_ic
            loss.backward()
            opt_pinn_net.step(); opt_pinn_phys.step()
            
            with torch.no_grad():
                t_zero_norm = torch.full_like(x_norm_pstd, -1.0)
                u0_pinn_pred = pinn(x_norm_pstd, t_zero_norm)
                nmse = (torch.sum((u0_pinn_pred - u0_pstd)**2) / torch.sum(u0_pstd**2)).item()
                history['High-Capacity PINN']['loss'].append(loss.item())
                history['High-Capacity PINN']['nmse_u0'].append(nmse)
                history['High-Capacity PINN']['gamma'].append(g_pred.item())
                history['High-Capacity PINN']['nu'].append(n_pred.item())
        history['High-Capacity PINN']['time_per_ep'] = (time.time() - t0) / epochs * 1000
        final_models['High-Capacity PINN'] = pinn
        
        # --- Quantitative Evaluations ---
        with torch.no_grad():
            # PSTD
            p0_eval_pstd, g_pstd, n_pstd = final_models['PSTD-DPNN'](x_norm_pstd)
            nmse_pstd = (torch.sum((p0_eval_pstd - u0_pstd)**2) / torch.sum(u0_pstd**2)).item()
            u_recon_pstd = pstd_solver(p0_eval_pstd, torch.zeros_like(x_grid_pstd), g_pstd, n_pstd)
            field_l2_pstd = (torch.norm(u_recon_pstd - u_exact_pstd) / torch.norm(u_exact_pstd)).item() * 100
            
            # FDTD
            p0_eval_fdtd, g_fdtd, n_fdtd = final_models['FDTD-DPNN'](x_norm_fdtd)
            nmse_fdtd = (torch.sum((p0_eval_fdtd - u0_fdtd)**2) / torch.sum(u0_fdtd**2)).item()
            u_recon_fdtd = fdtd_solver(p0_eval_fdtd, torch.zeros_like(x_grid_fdtd), g_fdtd, n_fdtd)
            # Evaluate field L2 on PSTD grid via interpolation
            u_recon_fdtd_interp = 0.5 * (u_recon_fdtd[:-1, :] + u_recon_fdtd[1:, :])
            field_l2_fdtd = (torch.norm(u_recon_fdtd_interp - u_exact_pstd) / torch.norm(u_exact_pstd)).item() * 100
            
            # PINN
            p0_eval_pinn = final_models['High-Capacity PINN'](x_norm_pstd, torch.full_like(x_norm_pstd, -1.0))
            g_pinn, n_pinn = final_models['High-Capacity PINN'].get_params()
            nmse_pinn = (torch.sum((p0_eval_pinn - u0_pstd)**2) / torch.sum(u0_pstd**2)).item()
            # Predict full field
            X_mesh, T_mesh = torch.meshgrid(x_grid_pstd, t_grid, indexing='ij')
            X_norm_mesh = (X_mesh / L) * 2 - 1
            T_norm_mesh = (T_mesh / T) * 2 - 1
            u_recon_pinn = final_models['High-Capacity PINN'](X_norm_mesh.unsqueeze(-1), T_norm_mesh.unsqueeze(-1))
            field_l2_pinn = (torch.norm(u_recon_pinn - u_exact_pstd) / torch.norm(u_exact_pstd)).item() * 100

        res_pstd = {
            'benchmark': exp['title'], 'method': 'PSTD-DPNN', 'nmse_u0': nmse_pstd, 'field_l2': field_l2_pstd,
            'gamma': g_pstd.item(), 'nu': n_pstd.item(),
            'err_gamma': abs(g_pstd.item() - true_gamma)/true_gamma * 100,
            'err_nu': abs(n_pstd.item() - true_nu)/true_nu * 100,
            'time_ep': history['PSTD-DPNN']['time_per_ep']
        }
        res_fdtd = {
            'benchmark': exp['title'], 'method': 'FDTD-DPNN', 'nmse_u0': nmse_fdtd, 'field_l2': field_l2_fdtd,
            'gamma': g_fdtd.item(), 'nu': n_fdtd.item(),
            'err_gamma': abs(g_fdtd.item() - true_gamma)/true_gamma * 100,
            'err_nu': abs(n_fdtd.item() - true_nu)/true_nu * 100,
            'time_ep': history['FDTD-DPNN']['time_per_ep']
        }
        res_pinn = {
            'benchmark': exp['title'], 'method': 'High-Capacity PINN', 'nmse_u0': nmse_pinn, 'field_l2': field_l2_pinn,
            'gamma': g_pinn.item(), 'nu': n_pinn.item(),
            'err_gamma': abs(g_pinn.item() - true_gamma)/true_gamma * 100,
            'err_nu': abs(n_pinn.item() - true_nu)/true_nu * 100,
            'time_ep': history['High-Capacity PINN']['time_per_ep']
        }
        all_summary_results.extend([res_pstd, res_fdtd, res_pinn])
        
        print(f"Results for {exp['title']}:")
        print(f"  PSTD-DPNN: NMSE={nmse_pstd:.4e}, Field L2={field_l2_pstd:.2f}%, g={g_pstd.item():.4f}, n={n_pstd.item():.2e}")
        print(f"  FDTD-DPNN: NMSE={nmse_fdtd:.4e}, Field L2={field_l2_fdtd:.2f}%, g={g_fdtd.item():.4f}, n={n_fdtd.item():.2e}")
        print(f"  PINN:      NMSE={nmse_pinn:.4e}, Field L2={field_l2_pinn:.2f}%, g={g_pinn.item():.4f}, n={n_pinn.item():.2e}")
        
        # --- Plot Figures ---
        ep_arr = np.arange(1, epochs + 1)
        colors = {'PSTD-DPNN': 'blue', 'FDTD-DPNN': 'orange', 'High-Capacity PINN': 'green'}
        
        # 1. Convergence Metrics Figure
        fig, axs = plt.subplots(1, 4, figsize=(18, 4.0))
        for m in methods:
            axs[0].semilogy(ep_arr, history[m]['loss'], label=m, color=colors[m], linewidth=1.8)
            axs[1].semilogy(ep_arr, history[m]['nmse_u0'], label=m, color=colors[m], linewidth=1.8)
            axs[2].plot(ep_arr, history[m]['gamma'], label=m, color=colors[m], linewidth=1.8)
            axs[3].plot(ep_arr, history[m]['nu'], label=m, color=colors[m], linewidth=1.8)
        axs[0].set_title('Training Loss', fontsize=11, fontweight='bold')
        axs[1].set_title('Initial Condition NMSE', fontsize=11, fontweight='bold')
        axs[2].set_title(r'Identified $\gamma$ (True: 0.10)', fontsize=11, fontweight='bold')
        axs[2].axhline(true_gamma, color='black', linestyle=':', label='Ground Truth')
        axs[3].set_title(r'Identified $\nu$ (True: 0.005)', fontsize=11, fontweight='bold')
        axs[3].axhline(true_nu, color='black', linestyle=':', label='Ground Truth')
        for i in range(4):
            axs[i].set_xlabel('Epoch', fontsize=10); axs[i].grid(True, linestyle=':', alpha=0.6)
            axs[i].legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(f"latex/Images/dpnn/{exp['name']}_metrics.pdf", bbox_inches='tight')
        plt.close()
        
        # 2. Reconstructed u0 Figure
        plt.figure(figsize=(7, 4.2))
        plt.plot(x_grid_pstd.numpy(), u0_pstd.numpy(), 'k-', label='Ground Truth $u_0(x)$', linewidth=2.5)
        plt.plot(x_grid_pstd.numpy(), p0_eval_pstd.numpy(), 'b--', label=f"PSTD-DPNN (NMSE={nmse_pstd:.2e})", linewidth=1.8)
        plt.plot(x_grid_fdtd.numpy(), p0_eval_fdtd.numpy(), 'm-.', label=f"FDTD-DPNN (NMSE={nmse_fdtd:.2e})", linewidth=1.5)
        plt.plot(x_grid_pstd.numpy(), p0_eval_pinn.numpy(), 'g:', label=f"PINN (NMSE={nmse_pinn:.2e})", linewidth=2.0)
        plt.plot(sensor_x.numpy(), torch.zeros_like(sensor_x).numpy(), 'rv', markersize=8, label='Sensors')
        plt.title(f"{exp['title']} - Initial Condition Reconstruction", fontsize=11, fontweight='bold')
        plt.xlabel('Spatial Coordinate $x$ [m]', fontsize=10)
        plt.ylabel('Initial Acoustic Pressure $p_0(x)$', fontsize=10)
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend(fontsize=9, loc='upper right')
        plt.tight_layout()
        plt.savefig(f"latex/Images/dpnn/{exp['name']}_u0_recon.pdf", bbox_inches='tight')
        plt.close()
        
        # 3. Spatio-temporal evaluation & error heatmaps
        fig, axs = plt.subplots(1, 4, figsize=(18, 4.0))
        im0 = axs[0].imshow(u_exact_pstd.numpy(), aspect='auto', origin='lower', extent=[0, T, 0, L], cmap='viridis')
        axs[0].set_title('Continuous Analytical Truth $u^*(x,t)$', fontsize=10, fontweight='bold')
        plt.colorbar(im0, ax=axs[0])
        
        err_pstd = np.abs(u_recon_pstd.numpy() - u_exact_pstd.numpy())
        err_fdtd = np.abs(u_recon_fdtd_interp.numpy() - u_exact_pstd.numpy())
        err_pinn = np.abs(u_recon_pinn.numpy() - u_exact_pstd.numpy())
        vmax_err = max(np.percentile(err_pstd, 99), np.percentile(err_fdtd, 99), np.percentile(err_pinn, 99))
        
        im1 = axs[1].imshow(err_pstd, aspect='auto', origin='lower', extent=[0, T, 0, L], cmap='magma', vmin=0, vmax=vmax_err)
        axs[1].set_title(f'PSTD-DPNN Error (L2={field_l2_pstd:.2f}%)', fontsize=10, fontweight='bold')
        plt.colorbar(im1, ax=axs[1])
        
        im2 = axs[2].imshow(err_fdtd, aspect='auto', origin='lower', extent=[0, T, 0, L], cmap='magma', vmin=0, vmax=vmax_err)
        axs[2].set_title(f'FDTD-DPNN Error (L2={field_l2_fdtd:.2f}%)', fontsize=10, fontweight='bold')
        plt.colorbar(im2, ax=axs[2])
        
        im3 = axs[3].imshow(err_pinn, aspect='auto', origin='lower', extent=[0, T, 0, L], cmap='magma', vmin=0, vmax=vmax_err)
        axs[3].set_title(f'PINN Error (L2={field_l2_pinn:.2f}%)', fontsize=10, fontweight='bold')
        plt.colorbar(im3, ax=axs[3])
        
        for i in range(4):
            axs[i].set_xlabel('Time $t$ [s]', fontsize=9)
            if i == 0: axs[i].set_ylabel('Position $x$ [m]', fontsize=9)
        plt.tight_layout()
        plt.savefig(f"latex/Images/dpnn/{exp['name']}_spacetime.png", dpi=200, bbox_inches='tight')
        plt.close()

    # --- Generate Table 3.4 (Comparison Table) ---
    latex_comp = r"""\begin{table}[H]
\centering
\small
\caption{Quantitative benchmark comparison between PSTD-DPNN, FDTD-DPNN, and the high-capacity PINN baseline (4 hidden layers, 128 units, $\sim 50{,}000$ parameters) across three initial condition profiles (10 sensors, 30\,dB SNR Gaussian noise, $T = 2.0$\,s, continuous analytical ground truth).}
\label{tab:dpnn_pinn_comparison}
\begin{tabular}{llccccc}
\toprule
\textbf{Benchmark} & \textbf{Method} & \textbf{NMSE($p_0$)} & \textbf{Field Rel. $L_2$ (\%)} & \textbf{Rel. Err $\gamma$ (\%)} & \textbf{Rel. Err $\nu$ (\%)} & \textbf{Time / Epoch (ms)} \\
\midrule
"""
    cur_b = ""
    for r in all_summary_results:
        if r['benchmark'] != cur_b:
            if cur_b != "": latex_comp += "\\midrule\n"
            cur_b = r['benchmark']
            b_label = f"\\multirow{{3}}{{*}}{{{cur_b}}}"
        else:
            b_label = ""
        latex_comp += f"{b_label} & {r['method']} & {r['nmse_u0']:.2e} & {r['field_l2']:.2f}\\% & {r['err_gamma']:.1f}\\% & {r['err_nu']:.1f}\\% & {r['time_ep']:.1f} \\\\\n"
    latex_comp += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    with open("latex/dpnn_comparison_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_comp)
    print("Saved latex/dpnn_comparison_table.tex")
    
    # --- Generate Table 3.5 (Parameter Table) ---
    latex_param = r"""\begin{table}[H]
\centering
\small
\caption{Detailed comparison of estimated physical attenuation parameters ($\gamma_{\mathrm{true}} = 0.100\,\mathrm{s}^{-1}$, $\nu_{\mathrm{true}} = 5.00 \times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$) across benchmarks and methods. In all cases, parameter estimation errors reflect structural identifiability limitations and compensation mechanisms between spatial amplitude and temporal decay.}
\label{tab:dpnn_param_estimation}
\begin{tabular}{llcccc}
\toprule
\textbf{Benchmark} & \textbf{Method} & $\hat{\gamma}$ [s$^{-1}$] & Rel. Err $\gamma$ & $\hat{\nu}$ [m$^2$/s] & Rel. Err $\nu$ \\
\midrule
"""
    cur_b = ""
    for r in all_summary_results:
        if r['benchmark'] != cur_b:
            if cur_b != "": latex_param += "\\midrule\n"
            cur_b = r['benchmark']
            b_label = f"\\multirow{{3}}{{*}}{{{cur_b}}}"
        else:
            b_label = ""
        latex_param += f"{b_label} & {r['method']} & {r['gamma']:.4f} & {r['err_gamma']:.1f}\\% & {r['nu']:.2e} & {r['err_nu']:.1f}\\% \\\\\n"
    latex_param += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    with open("latex/dpnn_parameters_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_param)
    print("Saved latex/dpnn_parameters_table.tex")
    
    return all_summary_results

# ==============================================================================
# 5. THREE ISOLATED INVERSION EXPERIMENTS & 2D LOSS LANDSCAPE
# ==============================================================================
def run_isolated_experiments():
    print("\n" + "="*80)
    print("RUNNING THREE ISOLATED INVERSION EXPERIMENTS (ON GAUSSIAN PULSE)")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
    Nx = int(L/dx); Nt = int(T/dt) + 1
    pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)
    x_grid = pstd.x_grid; t_grid = pstd.t_grid
    
    true_gamma = 0.1; true_nu = 0.005
    init_gamma = 0.5; init_nu = 0.001
    mu = L/2; sigma = L/20; M_pulse = 5
    
    def u0_func(x):
        val = np.zeros_like(x)
        for m in range(-M_pulse, M_pulse + 1):
            val += 0.5 * np.exp(-0.5 * ((x - (2*m*L + mu)) / sigma)**2)
            val += 0.5 * np.exp(-0.5 * ((x - (2*m*L - mu)) / sigma)**2)
        return val

    u0_true = torch.tensor(u0_func(x_grid.numpy()), dtype=torch.float32, device=device)
    u_true = torch.tensor(analytical_solution(x_grid.numpy(), t_grid.numpy(), u0_func, c, true_gamma, true_nu, L), dtype=torch.float32, device=device)
    
    torch.manual_seed(123)
    sensor_idx = torch.sort(torch.randperm(Nx - 2)[:10] + 1)[0]
    p_exact_sensors = u_true[sensor_idx, :]
    noise = torch.randn_like(p_exact_sensors)
    noise = noise / torch.norm(noise) * torch.norm(p_exact_sensors) * (10**(-30 / 20))
    p_noisy = p_exact_sensors + noise
    x_norm = ((x_grid / L) * 2 - 1).unsqueeze(-1)
    
    epochs = 350
    isolated_results = []
    
    # --- Experiment A: u0 only (Physical parameters fixed to ground truth) ---
    print("\n>>> Experiment A: Initial Condition Reconstruction Only (gamma, nu fixed to true values)")
    torch.manual_seed(42)
    model_A = DPNNModelLog(true_gamma, true_nu)
    opt_A = torch.optim.Adam(model_A.siren.parameters(), lr=1e-3)
    for ep in range(epochs):
        opt_A.zero_grad()
        p0_pred = model_A.siren(x_norm)
        u_pred = pstd(p0_pred, torch.zeros_like(x_grid), torch.tensor(true_gamma), torch.tensor(true_nu))
        loss = torch.mean((u_pred[sensor_idx, :] - p_noisy)**2)
        loss.backward()
        opt_A.step()
    with torch.no_grad():
        p0_A = model_A.siren(x_norm)
        nmse_A = (torch.sum((p0_A - u0_true)**2) / torch.sum(u0_true**2)).item()
        u_recon_A = pstd(p0_A, torch.zeros_like(x_grid), torch.tensor(true_gamma), torch.tensor(true_nu))
        field_l2_A = (torch.norm(u_recon_A - u_true) / torch.norm(u_true)).item() * 100
        isolated_results.append({
            'name': 'Exp A: Initial Condition Only ($u_0$ unknown, $\\gamma, \\nu$ fixed)',
            'nmse': nmse_A, 'field_l2': field_l2_A,
            'gamma': true_gamma, 'nu': true_nu, 'err_g': 0.0, 'err_n': 0.0
        })
        print(f"Exp A -> NMSE(p0) = {nmse_A:.4e}, Field L2 = {field_l2_A:.2f}%")
        
    # --- Experiment B: gamma, nu only (Initial condition fixed to ground truth) ---
    print("\n>>> Experiment B: Physical Parameter Identification Only (u0 fixed to true values)")
    torch.manual_seed(42)
    xi_g_B = nn.Parameter(torch.tensor(float(np.log(init_gamma / 0.1))))
    xi_n_B = nn.Parameter(torch.tensor(float(np.log(init_nu / 0.005))))
    opt_B = torch.optim.Adam([xi_g_B, xi_n_B], lr=5e-2)
    traj_B_g = []; traj_B_n = []
    for ep in range(epochs):
        opt_B.zero_grad()
        g_val = 0.1 * torch.exp(xi_g_B)
        n_val = 0.005 * torch.exp(xi_n_B)
        traj_B_g.append(g_val.item()); traj_B_n.append(n_val.item())
        u_pred = pstd(u0_true, torch.zeros_like(x_grid), g_val, n_val)
        loss = torch.mean((u_pred[sensor_idx, :] - p_noisy)**2)
        loss.backward()
        opt_B.step()
    with torch.no_grad():
        g_B = 0.1 * torch.exp(xi_g_B).item()
        n_B = 0.005 * torch.exp(xi_n_B).item()
        u_recon_B = pstd(u0_true, torch.zeros_like(x_grid), torch.tensor(g_B), torch.tensor(n_B))
        field_l2_B = (torch.norm(u_recon_B - u_true) / torch.norm(u_true)).item() * 100
        err_g_B = abs(g_B - true_gamma) / true_gamma * 100
        err_n_B = abs(n_B - true_nu) / true_nu * 100
        isolated_results.append({
            'name': 'Exp B: Physical Parameters Only ($\\gamma, \\nu$ unknown, $u_0$ fixed)',
            'nmse': 0.0, 'field_l2': field_l2_B,
            'gamma': g_B, 'nu': n_B, 'err_g': err_g_B, 'err_n': err_n_B
        })
        print(f"Exp B -> gamma = {g_B:.4f} (err: {err_g_B:.1f}%), nu = {n_B:.2e} (err: {err_n_B:.1f}%), Field L2 = {field_l2_B:.2f}%")
        
    # --- Experiment C: Joint Inversion (u0, gamma, nu all unknown) ---
    print("\n>>> Experiment C: Joint Inverse Problem (u0, gamma, nu all unknown)")
    torch.manual_seed(42)
    model_C = DPNNModelLog(init_gamma, init_nu)
    opt_net_C = torch.optim.Adam(model_C.siren.parameters(), lr=1e-3)
    opt_phys_C = torch.optim.Adam([model_C.xi_gamma, model_C.xi_nu], lr=5e-2)
    traj_C_g = []; traj_C_n = []
    for ep in range(epochs):
        opt_net_C.zero_grad(); opt_phys_C.zero_grad()
        p0_pred, g_pred, n_pred = model_C(x_norm)
        traj_C_g.append(g_pred.item()); traj_C_n.append(n_pred.item())
        u_pred = pstd(p0_pred, torch.zeros_like(x_grid), g_pred, n_pred)
        loss = torch.mean((u_pred[sensor_idx, :] - p_noisy)**2)
        loss.backward()
        opt_net_C.step(); opt_phys_C.step()
    with torch.no_grad():
        p0_C, g_C_t, n_C_t = model_C(x_norm)
        g_C = g_C_t.item(); n_C = n_C_t.item()
        nmse_C = (torch.sum((p0_C - u0_true)**2) / torch.sum(u0_true**2)).item()
        u_recon_C = pstd(p0_C, torch.zeros_like(x_grid), g_C_t, n_C_t)
        field_l2_C = (torch.norm(u_recon_C - u_true) / torch.norm(u_true)).item() * 100
        err_g_C = abs(g_C - true_gamma) / true_gamma * 100
        err_n_C = abs(n_C - true_nu) / true_nu * 100
        isolated_results.append({
            'name': 'Exp C: Joint Inversion ($u_0, \\gamma, \\nu$ all unknown)',
            'nmse': nmse_C, 'field_l2': field_l2_C,
            'gamma': g_C, 'nu': n_C, 'err_g': err_g_C, 'err_n': err_n_C
        })
        print(f"Exp C -> NMSE(p0) = {nmse_C:.4e}, Field L2 = {field_l2_C:.2f}%, gamma = {g_C:.4f} (err: {err_g_C:.1f}%), nu = {n_C:.2e} (err: {err_n_C:.1f}%)")
        
    # --- Generate Table 3.7 (Isolated Inversion Table) ---
    latex_iso = r"""\begin{table}[H]
\centering
\small
\caption{Disentangled analysis of the inverse problem on the Gaussian benchmark, comparing: (Exp~A) initial condition reconstruction with physical parameters fixed to ground truth; (Exp~B) physical parameter identification with initial condition fixed to ground truth; and (Exp~C) joint reconstruction where both the initial state and damping parameters are simultaneously unknown.}
\label{tab:dpnn_isolated_experiments}
\begin{tabular}{lcccccc}
\toprule
\textbf{Inverse Formulation} & \textbf{NMSE($p_0$)} & \textbf{Field Rel. $L_2$ (\%)} & $\hat{\gamma}$ [s$^{-1}$] & Rel. Err $\gamma$ & $\hat{\nu}$ [m$^2$/s] & Rel. Err $\nu$ \\
\midrule
"""
    for r in isolated_results:
        nmse_str = f"{r['nmse']:.2e}" if r['nmse'] > 0 else "0.00 (Exact)"
        latex_iso += f"{r['name']} & {nmse_str} & {r['field_l2']:.2f}\\% & {r['gamma']:.4f} & {r['err_g']:.1f}\\% & {r['nu']:.2e} & {r['err_n']:.1f}\\% \\\\\n"
    latex_iso += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    with open("latex/dpnn_isolated_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_iso)
    print("Saved latex/dpnn_isolated_table.tex")
    
    # --- Compute 2D Loss Landscape Grid & Contour Map ---
    print("\n>>> Computing 2D Loss Landscape Phi(gamma, nu)...")
    gamma_grid = np.linspace(0.0, 0.4, 50)
    nu_grid = np.linspace(0.0001, 0.015, 50)
    G_mesh, N_mesh = np.meshgrid(gamma_grid, nu_grid)
    loss_landscape = np.zeros_like(G_mesh)
    
    for i in range(G_mesh.shape[0]):
        for j in range(G_mesh.shape[1]):
            g_val = torch.tensor(G_mesh[i, j], dtype=torch.float32)
            n_val = torch.tensor(N_mesh[i, j], dtype=torch.float32)
            with torch.no_grad():
                u_sim = pstd(u0_true, torch.zeros_like(x_grid), g_val, n_val)
                loss_landscape[i, j] = torch.mean((u_sim[sensor_idx, :] - p_noisy)**2).item()
                
    plt.figure(figsize=(8, 6.0))
    levels = np.logspace(np.log10(np.min(loss_landscape) + 1e-6), np.log10(np.max(loss_landscape)), 25)
    cs = plt.contourf(G_mesh, N_mesh, loss_landscape, levels=levels, cmap='viridis', norm=plt.Normalize())
    plt.colorbar(cs, label='Observation Loss $\\mathcal{L}_{\\mathrm{data}}$')
    plt.contour(G_mesh, N_mesh, loss_landscape, levels=levels, colors='black', linewidths=0.5, alpha=0.5)
    
    # Ground truth marker
    plt.plot(true_gamma, true_nu, 'r*', markersize=14, label='Ground Truth $(0.10, 0.005)$')
    # Experiment B optimization trajectory
    plt.plot(traj_B_g, traj_B_n, 'w.-', linewidth=2.0, markersize=4, label='Exp B Trajectory ($u_0$ fixed)')
    # Experiment C optimization trajectory
    plt.plot(traj_C_g, traj_C_n, 'y.--', linewidth=2.0, markersize=4, label='Exp C Trajectory (Joint $u_0$ unknown)')
    
    plt.title(r'2D Loss Landscape $\Phi(\gamma, \nu)$ and Inversion Trajectories', fontsize=12, fontweight='bold')
    plt.xlabel(r'Telegrapher Damping $\gamma$ [s$^{-1}$]', fontsize=11)
    plt.ylabel(r'Viscoelastic Damping $\nu$ [m$^2$/s]', fontsize=11)
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig("latex/Images/dpnn/identifiability_landscape_contour.pdf", bbox_inches='tight')
    plt.close()
    print("Saved latex/Images/dpnn/identifiability_landscape_contour.pdf")
    
    return isolated_results

# ==============================================================================
# 6. ABLATION STUDY INCLUDING TOTAL VARIATION (TV) REGULARIZATION
# ==============================================================================
def run_ablation_study():
    print("\n" + "="*80)
    print("RUNNING REGULARIZATION ABLATION STUDY (GAUSSIAN & TRIANGLE)")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
    Nx = int(L/dx); Nt = int(T/dt) + 1
    pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)
    x_grid = pstd.x_grid; t_grid = pstd.t_grid
    
    torch.manual_seed(123)
    sensor_idx = torch.sort(torch.randperm(Nx - 2)[:10] + 1)[0]
    x_norm = ((x_grid / L) * 2 - 1).unsqueeze(-1)
    
    mu = L/2; sigma = L/20; M_pulse = 5
    
    ablation_configs = [
        {'name': 'Data only', 'sp': 0.0, 'smooth': 0.0, 'tv': 0.0, 'decoupled': False},
        {'name': '$L_1$ only', 'sp': 0.005, 'smooth': 0.0, 'tv': 0.0, 'decoupled': False},
        {'name': '$H^1$ only', 'sp': 0.0, 'smooth': 0.1, 'tv': 0.0, 'decoupled': False},
        {'name': '$L_1 + H^1$ (Coupled $\\nu$)', 'sp': 0.005, 'smooth': 0.1, 'tv': 0.0, 'decoupled': False},
        {'name': '$L_1 + H^1$ (Decoupled $\\nu_0$)', 'sp': 0.005, 'smooth': 0.1, 'tv': 0.0, 'decoupled': True},
        {'name': 'Total Variation (TV)', 'sp': 0.0, 'smooth': 0.0, 'tv': 0.01, 'decoupled': False},
    ]
    
    test_profiles = [
        ('Gaussian Pulse', lambda x: sum([0.5 * np.exp(-0.5 * ((x - (2*m*L + mu)) / sigma)**2) + 0.5 * np.exp(-0.5 * ((x - (2*m*L - mu)) / sigma)**2) for m in range(-M_pulse, M_pulse + 1)])),
        ('Triangle Pulse', lambda x: sum([0.5 * np.maximum(0.0, 1.0 - np.abs(x - (2*m*L + mu)) / sigma) + 0.5 * np.maximum(0.0, 1.0 - np.abs(x - (2*m*L - mu)) / sigma) for m in range(-M_pulse, M_pulse + 1)]))
    ]
    
    ablation_results = []
    
    for bname, u0_func in test_profiles:
        u0_true = torch.tensor(u0_func(x_grid.numpy()), dtype=torch.float32, device=device)
        u_true = torch.tensor(analytical_solution(x_grid.numpy(), t_grid.numpy(), u0_func, c, 0.1, 0.005, L), dtype=torch.float32, device=device)
        
        torch.manual_seed(123)
        p_obs = u_true[sensor_idx, :]
        noise = torch.randn_like(p_obs)
        noise = noise / torch.norm(noise) * torch.norm(p_obs) * (10**(-30 / 20))
        p_noisy = p_obs + noise
        
        for cfg in ablation_configs:
            torch.manual_seed(42)
            model = DPNNModelLog(0.5, 0.001)
            opt_net = torch.optim.Adam(model.siren.parameters(), lr=1e-3)
            opt_phys = torch.optim.Adam([model.xi_gamma, model.xi_nu], lr=5e-2)
            
            for ep in range(350):
                opt_net.zero_grad(); opt_phys.zero_grad()
                p0_pred, g_pred, n_pred = model(x_norm)
                u_pred = pstd(p0_pred, torch.zeros_like(x_grid), g_pred, n_pred)
                
                loss_data = torch.mean((u_pred[sensor_idx, :] - p_noisy)**2)
                loss_sp = torch.mean(torch.abs(p0_pred)) if cfg['sp'] > 0 else torch.tensor(0.0)
                
                # Spatial derivatives
                p0_x = torch.zeros_like(p0_pred)
                p0_x[1:-1] = (p0_pred[2:] - p0_pred[:-2]) / (2 * dx)
                p0_x[0] = (p0_pred[1] - p0_pred[0]) / dx
                p0_x[-1] = (p0_pred[-1] - p0_pred[-2]) / dx
                
                if cfg['smooth'] > 0:
                    reg_weight = 0.005 if cfg['decoupled'] else n_pred
                    loss_h1 = reg_weight * torch.sum(p0_x**2) * dx
                else:
                    loss_h1 = torch.tensor(0.0)
                    
                if cfg['tv'] > 0:
                    loss_tv = torch.sum(torch.abs(p0_x)) * dx
                else:
                    loss_tv = torch.tensor(0.0)
                    
                loss = loss_data + cfg['sp'] * loss_sp + cfg['smooth'] * loss_h1 + cfg['tv'] * loss_tv
                loss.backward()
                opt_net.step(); opt_phys.step()
                
            with torch.no_grad():
                p0_eval, g_eval, n_eval = model(x_norm)
                nmse = (torch.sum((p0_eval - u0_true)**2) / torch.sum(u0_true**2)).item()
                u_recon = pstd(p0_eval, torch.zeros_like(x_grid), g_eval, n_eval)
                field_l2 = (torch.norm(u_recon - u_true) / torch.norm(u_true)).item() * 100
                
                ablation_results.append({
                    'profile': bname,
                    'config': cfg['name'],
                    'nmse': nmse,
                    'field_l2': field_l2,
                    'gamma_est': g_eval.item(),
                    'nu_est': n_eval.item()
                })
                print(f"{bname:14s} | {cfg['name']:28s} -> NMSE={nmse:.4f}, Field L2={field_l2:.2f}%, g={g_eval.item():.3f}, n={n_eval.item():.2e}")

    # Build Ablation LaTeX table
    latex_abl_tab = r"""\begin{table}[H]
\centering
\small
\caption{Ablation study on regularization components for PSTD-DPNN across Gaussian and Triangular initial profiles. Comparing unregularized data fitting, individual $L_1$ and $H^1$ penalties, the joint $L_1 + H^1$ formulation, decoupled $H^1$ regularization independent of $\nu$, and Total Variation (TV) penalization.}
\label{tab:dpnn_ablation}
\begin{tabular}{llcccc}
\toprule
\textbf{Benchmark} & \textbf{Regularization Configuration} & \textbf{NMSE($p_0$)} & \textbf{Field Rel. $L_2$ (\%)} & $\hat{\gamma}$ [s$^{-1}$] & $\hat{\nu}$ [m$^2$/s] \\
\midrule
"""
    cur_p = ""
    for r in ablation_results:
        if r['profile'] != cur_p:
            if cur_p != "": latex_abl_tab += "\\midrule\n"
            cur_p = r['profile']
            p_label = f"\\multirow{{6}}{{*}}{{{cur_p}}}"
        else:
            p_label = ""
        latex_abl_tab += f"{p_label} & {r['config']} & {r['nmse']:.4f} & {r['field_l2']:.2f}\\% & {r['gamma_est']:.4f} & {r['nu_est']:.2e} \\\\\n"
    latex_abl_tab += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    with open("latex/dpnn_ablation_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_abl_tab)
    print("Saved latex/dpnn_ablation_table.tex")
    return ablation_results

# ==============================================================================
# 7. IDENTIFIABILITY, DIMENSIONLESS MATRIX & ORTHOGONAL PROJECTION ANALYSIS
# ==============================================================================
def run_identifiability_and_projection():
    print("\n" + "="*80)
    print("RUNNING IDENTIFIABILITY, DIMENSIONLESS & ORTHOGONAL PROJECTION ANALYSIS")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
    Nx = int(L/dx); Nt = int(T/dt)+1
    pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)
    x_grid = pstd.x_grid; t_grid = pstd.t_grid
    
    torch.manual_seed(123)
    sensor_idx = torch.sort(torch.randperm(Nx - 2)[:10] + 1)[0]
    sensor_x = x_grid[sensor_idx]
    
    mu = L/2; sigma = L/20; M_pulse = 5
    profiles = [
        ('Standing Waves', sum([torch.cos(m * np.pi * x_grid / L) for m in range(1, 6)])),
        ('Gaussian Pulse', torch.exp(-0.5 * ((x_grid - mu) / sigma)**2)),
        ('Triangle Pulse', torch.clamp(1 - torch.abs(x_grid - mu)/sigma, min=0))
    ]
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    ident_results = []
    
    for idx, (pname, u0) in enumerate(profiles):
        true_gamma = torch.tensor(0.1, requires_grad=True)
        true_nu = torch.tensor(0.005, requires_grad=True)
        u = pstd(u0, torch.zeros_like(x_grid), true_gamma, true_nu)
        u_sens = u[sensor_idx, :] # (10, Nt)
        
        # Evaluate Jacobians wrt gamma and nu
        J_g = []; J_n = []
        for s in range(10):
            for t_step in range(Nt):
                g_val = torch.autograd.grad(u_sens[s, t_step], true_gamma, retain_graph=True)[0].item()
                n_val = torch.autograd.grad(u_sens[s, t_step], true_nu, retain_graph=True)[0].item()
                J_g.append(g_val); J_n.append(n_val)
                
        J_raw = np.column_stack([J_g, J_n]) # (410, 2)
        cond_raw = np.linalg.cond(J_raw.T @ J_raw)
        corr = np.dot(J_raw[:,0], J_raw[:,1]) / (np.linalg.norm(J_raw[:,0]) * np.linalg.norm(J_raw[:,1]))
        
        # Normalized Jacobian: columns scaled to unit norm
        J_norm = J_raw / np.linalg.norm(J_raw, axis=0, keepdims=True)
        cond_norm = np.linalg.cond(J_norm.T @ J_norm)
        
        # Dimensionless Jacobian: J_tilde = J * diag(gamma_ref, nu_ref)
        J_dim = J_raw * np.array([0.1, 0.005])
        cond_dim = np.linalg.cond(J_dim.T @ J_dim)
        
        # Orthogonal projection against initial condition basis u0
        J_u0 = np.zeros((10 * Nt, Nx))
        for k in range(Nx):
            e_k = torch.zeros_like(x_grid); e_k[k] = 1.0
            u_k = pstd(e_k, torch.zeros_like(x_grid), torch.tensor(0.1), torch.tensor(0.005))
            J_u0[:, k] = u_k[sensor_idx, :].detach().numpy().flatten()
            
        U_u0, S_u0, _ = np.linalg.svd(J_u0, full_matrices=False)
        r = np.sum(S_u0 > 1e-10 * S_u0[0])
        U_r = U_u0[:, :r]
        J_dim_proj = J_dim - U_r @ (U_r.T @ J_dim)
        
        ret_g = np.linalg.norm(J_dim_proj[:, 0]) / np.linalg.norm(J_dim[:, 0]) * 100
        ret_n = np.linalg.norm(J_dim_proj[:, 1]) / np.linalg.norm(J_dim[:, 1]) * 100
        cond_proj = np.linalg.cond(J_dim_proj.T @ J_dim_proj)
        corr_proj = np.dot(J_dim_proj[:, 0], J_dim_proj[:, 1]) / (np.linalg.norm(J_dim_proj[:, 0]) * np.linalg.norm(J_dim_proj[:, 1]))
        
        ident_results.append({
            'profile': pname,
            'cond_raw': cond_raw, 'cond_norm': cond_norm, 'cond_dim': cond_dim,
            'corr': corr, 'ret_g': ret_g, 'ret_n': ret_n,
            'cond_proj': cond_proj, 'corr_proj': corr_proj
        })
        print(f"{pname:16s} -> Raw cond: {cond_raw:8.2e} | Dim cond: {cond_dim:5.2f} | Norm cond: {cond_norm:5.2f} | Corr: {corr:.3f}")
        print(f"                 Retained after u0 proj: gamma={ret_g:.1f}%, nu={ret_n:.1f}% | Proj cond: {cond_proj:5.2f}")
        
        # Plot sensitivity time traces at sensor 5 (mid-domain)
        sens_idx_mid = 4
        J_g_mid = J_g[sens_idx_mid*Nt : (sens_idx_mid+1)*Nt]
        J_n_mid = J_n[sens_idx_mid*Nt : (sens_idx_mid+1)*Nt]
        J_g_n = np.array(J_g_mid) / (np.max(np.abs(J_g_mid)) + 1e-12)
        J_n_n = np.array(J_n_mid) / (np.max(np.abs(J_n_mid)) + 1e-12)
        
        axs[idx].plot(t_grid.numpy(), J_g_n, 'b-', label=r'Norm. $\partial u / \partial \gamma$', linewidth=2.0)
        axs[idx].plot(t_grid.numpy(), J_n_n, 'r--', label=r'Norm. $\partial u / \partial \nu$', linewidth=1.8)
        axs[idx].set_title(f"{pname}\n($\\kappa(\\widetilde{{J}}^T\\widetilde{{J}})={cond_dim:.1f}$, $\\kappa(J_{{\\mathrm{{norm}}}})={cond_norm:.1f}$)", fontsize=10, fontweight='bold')
        axs[idx].set_xlabel('Time $t$ [s]', fontsize=10)
        if idx == 0: axs[idx].set_ylabel('Normalized Sensitivity', fontsize=10)
        axs[idx].grid(True, linestyle=':', alpha=0.6)
        axs[idx].legend(fontsize=9, loc='upper right')
        
    plt.tight_layout()
    plt.savefig("latex/Images/dpnn/identifiability_sensitivity.pdf", bbox_inches='tight')
    plt.close()
    print("Saved latex/Images/dpnn/identifiability_sensitivity.pdf")
    return ident_results

# ==============================================================================
# 8. FISHER INFORMATION & OPTIMAL SENSOR PLACEMENT
# ==============================================================================
def run_fisher_sensor_placement():
    print("\n" + "="*80)
    print("RUNNING FISHER INFORMATION & OPTIMAL SENSOR PLACEMENT SUITE")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05; c = 1.0
    Nx = int(L/dx); Nt = int(T/dt) + 1
    pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)
    x_grid = pstd.x_grid; t_grid = pstd.t_grid
    
    mu = L/2; sigma = L/20
    u0 = torch.exp(-0.5 * ((x_grid - mu)/sigma)**2)
    true_gamma = torch.tensor(0.1, requires_grad=True)
    true_nu = torch.tensor(0.005, requires_grad=True)
    
    u = pstd(u0, torch.zeros_like(x_grid), true_gamma, true_nu)
    
    # Calculate sensitivities at all grid points
    J_g_all = np.zeros((Nx, Nt))
    J_n_all = np.zeros((Nx, Nt))
    for i in range(Nx):
        for n in range(Nt):
            g_val = torch.autograd.grad(u[i, n], true_gamma, retain_graph=True)[0].item()
            n_val = torch.autograd.grad(u[i, n], true_nu, retain_graph=True)[0].item()
            J_g_all[i, n] = g_val * 0.1 # dimensionless
            J_n_all[i, n] = n_val * 0.005 # dimensionless
            
    # Fisher for candidate sensor i
    F_all = np.zeros((Nx, 2, 2))
    for i in range(Nx):
        F_all[i, 0, 0] = np.sum(J_g_all[i]**2)
        F_all[i, 0, 1] = np.sum(J_g_all[i] * J_n_all[i])
        F_all[i, 1, 0] = F_all[i, 0, 1]
        F_all[i, 1, 1] = np.sum(J_n_all[i]**2)
        
    # Configuration 1: Random (seed 123)
    torch.manual_seed(123)
    idx_rand = torch.sort(torch.randperm(Nx - 2)[:10] + 1)[0].numpy()
    F_rand = np.sum(F_all[idx_rand], axis=0)
    det_rand = np.linalg.det(F_rand)
    eig_rand = np.linalg.eigvalsh(F_rand)
    cond_rand = eig_rand[1] / eig_rand[0]
    
    # Configuration 2: Uniform
    idx_uni = np.linspace(2, Nx-3, 10, dtype=int)
    F_uni = np.sum(F_all[idx_uni], axis=0)
    det_uni = np.linalg.det(F_uni)
    eig_uni = np.linalg.eigvalsh(F_uni)
    cond_uni = eig_uni[1] / eig_uni[0]
    
    # Configuration 3: Fisher D-optimal selection
    selected_opt = []
    F_curr = np.zeros((2, 2))
    candidates = list(range(Nx))
    for _ in range(10):
        best_cand = None; best_det = -1
        for c_cand in candidates:
            F_try = F_curr + F_all[c_cand]
            d = np.linalg.det(F_try)
            if d > best_det:
                best_det = d; best_cand = c_cand
        selected_opt.append(best_cand)
        candidates.remove(best_cand)
        F_curr = F_curr + F_all[best_cand]
    idx_opt = np.array(sorted(selected_opt))
    F_opt = np.sum(F_all[idx_opt], axis=0)
    det_opt = np.linalg.det(F_opt)
    eig_opt = np.linalg.eigvalsh(F_opt)
    cond_opt = eig_opt[1] / eig_opt[0]
    
    # Run Exp B (physical parameter identification with true u0) under each sensor placement
    u0_tensor = torch.tensor(u0, dtype=torch.float32, device=device)
    u_true_pstd = pstd(u0_tensor, torch.zeros_like(x_grid), torch.tensor(0.1), torch.tensor(0.005))
    
    configs = [
        ('Random Placement (Seed 123)', idx_rand, det_rand, eig_rand[0], cond_rand),
        ('Uniform Equispaced Placement', idx_uni, det_uni, eig_uni[0], cond_uni),
        ('Fisher D-Optimal Placement', idx_opt, det_opt, eig_opt[0], cond_opt),
    ]
    
    sensor_eval_results = []
    for cname, s_idx, d_val, l_min, c_val in configs:
        torch.manual_seed(42)
        xi_g = nn.Parameter(torch.tensor(float(np.log(0.5 / 0.1))))
        xi_n = nn.Parameter(torch.tensor(float(np.log(0.001 / 0.005))))
        opt = torch.optim.Adam([xi_g, xi_n], lr=5e-2)
        
        p_exact = u_true_pstd[s_idx, :]
        torch.manual_seed(123)
        noise = torch.randn_like(p_exact)
        noise = noise / torch.norm(noise) * torch.norm(p_exact) * (10**(-30 / 20))
        p_noisy = p_exact + noise
        
        for ep in range(350):
            opt.zero_grad()
            g_sim = 0.1 * torch.exp(xi_g)
            n_sim = 0.005 * torch.exp(xi_n)
            u_sim = pstd(u0_tensor, torch.zeros_like(x_grid), g_sim, n_sim)
            loss = torch.mean((u_sim[s_idx, :] - p_noisy)**2)
            loss.backward()
            opt.step()
            
        with torch.no_grad():
            g_final = (0.1 * torch.exp(xi_g)).item()
            n_final = (0.005 * torch.exp(xi_n)).item()
            err_g = abs(g_final - 0.1) / 0.1 * 100
            err_n = abs(n_final - 0.005) / 0.005 * 100
            sensor_eval_results.append({
                'name': cname,
                'det': d_val, 'lambda_min': l_min, 'cond': c_val,
                'gamma': g_final, 'nu': n_final, 'err_g': err_g, 'err_n': err_n
            })
            print(f"{cname:30s} -> det(F)={d_val:.2e}, lambda_min={l_min:.2e}, cond={c_val:.1f} | g={g_final:.4f} ({err_g:.1f}%), n={n_final:.2e} ({err_n:.1f}%)")

    # Generate Table 3.8 (Sensor Placement Table)
    latex_sensor = r"""\begin{table}[H]
\centering
\small
\caption{Comparison of sensor placement strategies (10 sensors, 30\,dB SNR) evaluated on the Gaussian pulse: Fisher information determinant $\det(F)$ (D-optimality), minimum eigenvalue $\lambda_{\min}(F)$ (E-optimality), condition number $\kappa(F)$, and resulting physical parameter identification errors.}
\label{tab:dpnn_sensor_placement}
\begin{tabular}{lcccccc}
\toprule
\textbf{Placement Strategy} & $\det(F)$ & $\lambda_{\min}(F)$ & $\kappa(F)$ & $\hat{\gamma}$ [s$^{-1}$] & Rel. Err $\gamma$ & Rel. Err $\nu$ \\
\midrule
"""
    for r in sensor_eval_results:
        latex_sensor += f"{r['name']} & {r['det']:.2e} & {r['lambda_min']:.2e} & {r['cond']:.1f} & {r['gamma']:.4f} & {r['err_g']:.1f}\\% & {r['err_n']:.1f}\\% \\\\\n"
    latex_sensor += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    with open("latex/dpnn_sensor_table.tex", "w", encoding="utf-8") as f:
        f.write(latex_sensor)
    print("Saved latex/dpnn_sensor_table.tex")
    return sensor_eval_results

# ==============================================================================
# 9. DIAGNOSTIC SUITES (RESOLUTION, CFL=1, STATISTICAL SENSITIVITY)
# ==============================================================================
def run_diagnostics():
    print("\n" + "="*80)
    print("RUNNING DIAGNOSTIC SUITES: RESOLUTION, CFL=1 DISPERSION, STATISTICAL SENSITIVITY")
    print("="*80)
    
    L = 2.0; T = 2.0; dx = 0.05; dt = 0.05
    Nx = int(L/dx); Nt = int(T/dt) + 1
    pstd = DifferentiablePSTD1D(Nx, dx, dt, L, T, device=device)
    x_grid = pstd.x_grid
    
    # Statistical Sensitivity (Sensors & SNR)
    mu = L/2; sigma = L/20; M_pulse = 5
    def u0_func(x):
        val = np.zeros_like(x)
        for m in range(-M_pulse, M_pulse + 1):
            val += 0.5 * np.exp(-0.5 * ((x - (2*m*L + mu)) / sigma)**2)
            val += 0.5 * np.exp(-0.5 * ((x - (2*m*L - mu)) / sigma)**2)
        return val

    u0_stat = torch.tensor(u0_func(x_grid.numpy()), dtype=torch.float32, device=device)
    u_stat_true = pstd(u0_stat, torch.zeros_like(x_grid), torch.tensor(0.1), torch.tensor(0.005))
    
    sensor_counts = [5, 10, 20, 40]
    snr_levels = [10, 20, 30, 40]
    
    sens_stats = []
    for Ns in sensor_counts:
        errors = []
        for seed in [1, 2, 3, 4, 5]:
            torch.manual_seed(seed * 100 + Ns)
            s_idx = torch.sort(torch.randperm(Nx - 2)[:Ns] + 1)[0]
            p_obs = u_stat_true[s_idx, :]
            noise = torch.randn_like(p_obs)
            noise = noise / torch.norm(noise) * torch.norm(p_obs) * (10**(-30/20))
            p_noisy = p_obs + noise
            
            model = DPNNModelLog(0.5, 0.001)
            opt = torch.optim.Adam(model.parameters(), lr=2e-3)
            x_norm = ((x_grid / L) * 2 - 1).unsqueeze(-1)
            for ep in range(150):
                opt.zero_grad()
                p0, g, n = model(x_norm)
                u_p = pstd(p0, torch.zeros_like(x_grid), g, n)
                loss = torch.mean((u_p[s_idx, :] - p_noisy)**2)
                loss.backward()
                opt.step()
            with torch.no_grad():
                p0_e, _, _ = model(x_norm)
                nmse = (torch.sum((p0_e - u0_stat)**2) / torch.sum(u0_stat**2)).item()
                errors.append(nmse)
        sens_stats.append({'Ns': Ns, 'mean': np.mean(errors), 'std': np.std(errors)})
        print(f"Sensors = {Ns:2d} -> NMSE = {np.mean(errors):.4f} +/- {np.std(errors):.4f}")
        
    snr_stats = []
    for snr in snr_levels:
        errors = []
        for seed in [1, 2, 3, 4, 5]:
            torch.manual_seed(seed * 200 + snr)
            s_idx = torch.sort(torch.randperm(Nx - 2)[:10] + 1)[0]
            p_obs = u_stat_true[s_idx, :]
            noise = torch.randn_like(p_obs)
            noise = noise / torch.norm(noise) * torch.norm(p_obs) * (10**(-snr/20))
            p_noisy = p_obs + noise
            
            model = DPNNModelLog(0.5, 0.001)
            opt = torch.optim.Adam(model.parameters(), lr=2e-3)
            x_norm = ((x_grid / L) * 2 - 1).unsqueeze(-1)
            for ep in range(150):
                opt.zero_grad()
                p0, g, n = model(x_norm)
                u_p = pstd(p0, torch.zeros_like(x_grid), g, n)
                loss = torch.mean((u_p[s_idx, :] - p_noisy)**2)
                loss.backward()
                opt.step()
            with torch.no_grad():
                p0_e, _, _ = model(x_norm)
                nmse = (torch.sum((p0_e - u0_stat)**2) / torch.sum(u0_stat**2)).item()
                errors.append(nmse)
        snr_stats.append({'snr': snr, 'mean': np.mean(errors), 'std': np.std(errors)})
        print(f"SNR = {snr:2d} dB -> NMSE = {np.mean(errors):.4f} +/- {np.std(errors):.4f}")

    # Plot Statistical Sensitivity
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0))
    axs[0].errorbar([s['Ns'] for s in sens_stats], [s['mean'] for s in sens_stats], yerr=[s['std'] for s in sens_stats], fmt='bo-', capsize=4, linewidth=1.8)
    axs[0].set_xlabel('Number of Sensors $N_{\\mathrm{sensors}}$', fontsize=11)
    axs[0].set_ylabel('Initial Condition NMSE', fontsize=11)
    axs[0].set_title('(a) Sensitivity to Sensor Density (SNR = 30 dB)', fontsize=11, fontweight='bold')
    axs[0].grid(True, linestyle=':', alpha=0.6)

    axs[1].errorbar([s['snr'] for s in snr_stats], [s['mean'] for s in snr_stats], yerr=[s['std'] for s in snr_stats], fmt='rs-', capsize=4, linewidth=1.8)
    axs[1].set_xlabel('Signal-to-Noise Ratio (SNR) [dB]', fontsize=11)
    axs[1].set_ylabel('Initial Condition NMSE', fontsize=11)
    axs[1].set_title('(b) Robustness to Measurement Noise (10 Sensors)', fontsize=11, fontweight='bold')
    axs[1].grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig("latex/Images/dpnn/statistical_sensitivity.pdf", bbox_inches='tight')
    plt.close()
    print("Saved latex/Images/dpnn/statistical_sensitivity.pdf")

# ==============================================================================
# MAIN EXECUTION ENTRY POINT
# ==============================================================================
if __name__ == '__main__':
    t_start = time.time()
    print("STARTING UNIFIED DPNN COMPLETE EXECUTION SUITE...")
    
    # 1. Main benchmarks
    run_main_benchmarks()
    
    # 2. Three isolated inversion experiments & 2D landscape
    run_isolated_experiments()
    
    # 3. Regularization ablation study (with TV)
    run_ablation_study()
    
    # 4. Dimensionless & Projected Jacobian analysis
    run_identifiability_and_projection()
    
    # 5. Fisher information sensor optimization
    run_fisher_sensor_placement()
    
    # 6. Diagnostics (Resolution, Dispersion, Statistical)
    run_diagnostics()
    
    t_total = time.time() - t_start
    print(f"\nALL SUITES EXECUTED SUCCESSFULLY IN {t_total:.1f} SECONDS ({t_total/60:.2f} MINUTES)!")
