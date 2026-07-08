% test_G_matrix.m
clear all; close all; clc;
addpath(genpath('../src'));
addpath('WASAbi1D/utils');

Nx = 100;
L = 1;
dx = L/Nx;
c = 1;
gamma = 0;
nu = 0;

space_order = 2;
C_res = get_residue_matrix(Nx, space_order);

NxSub = Nx/2;
Lsub = L/2;
km = (0:NxSub-1)' * pi / Lsub;
lambda = -(km.^2);

state_size = 2*Nx;

% Construct A exactly as before
A = zeros(Nx, Nx);
for col = 1:Nx
    U_test = zeros(Nx, 1);
    U_test(col) = 1;
    
    UL = U_test(1:NxSub);
    UR = U_test(NxSub+1:Nx);
    
    uL = idct(UL);
    uR = idct(UR);
    
    u = [uL; uR];
    f = (c^2/dx^2) * (C_res * u);
    
    FL = dct(f(1:NxSub));
    FR = dct(f(NxSub+1:end));
    
    A(1:NxSub, col) = FL;
    A(NxSub+1:end, col) = FR;
end

CFL = 1.0;
dt = CFL * dx / c;

% Construct D exactly as before
[S_cells, Tm] = modal_propagators(lambda, dt, c, gamma, nu);
D = zeros(state_size, state_size);
for j = 1:NxSub
    S_mat = S_cells{j};
    D(j, j) = S_mat(1,1);
    D(j, Nx+j) = S_mat(1,2);
    D(Nx+j, j) = S_mat(2,1);
    D(Nx+j, Nx+j) = S_mat(2,2);
    
    D(NxSub+j, NxSub+j) = S_mat(1,1);
    D(NxSub+j, Nx+NxSub+j) = S_mat(1,2);
    D(Nx+NxSub+j, NxSub+j) = S_mat(2,1);
    D(Nx+NxSub+j, Nx+NxSub+j) = S_mat(2,2);
end

K = [eye(Nx),   zeros(Nx);
     (dt/2)*A,  eye(Nx)];

G = K * D * K;

% --- Random State Test ---
rng(42); % Reproducibility
Wn = randn(state_size, 1);
W_matrix = G * Wn;

% --- Exact Solver Steps ---
UL = Wn(1:NxSub);
UR = Wn(NxSub+1:Nx);
VL = Wn(Nx+1:Nx+NxSub);
VR = Wn(Nx+NxSub+1:end);

% First half-kick
uL_phys = idct(UL);
uR_phys = idct(UR);
u_phys = [uL_phys; uR_phys];
fCorr_global = (c^2/dx^2) * (C_res * u_phys);
fCorrL = fCorr_global(1:NxSub);
fCorrR = fCorr_global(NxSub+1:end);

FcorrL = dct(fCorrL);
FcorrR = dct(fCorrR);

VL = VL + (dt/2) * FcorrL;
VR = VR + (dt/2) * FcorrR;

% Drift (modal step)
for j = 1:NxSub
    zL = [UL(j); VL(j)];
    zL = S_cells{j}*zL + Tm{j}*[0; 0]; % FmL=0 during drift
    UL(j) = zL(1);
    VL(j) = zL(2);
    
    zR = [UR(j); VR(j)];
    zR = S_cells{j}*zR + Tm{j}*[0; 0];
    UR(j) = zR(1);
    VR(j) = zR(2);
end

% Second half-kick
uL_phys = idct(UL);
uR_phys = idct(UR);
u_phys = [uL_phys; uR_phys];
fCorr_global = (c^2/dx^2) * (C_res * u_phys);
fCorrL = fCorr_global(1:NxSub);
fCorrR = fCorr_global(NxSub+1:end);

FcorrL_new = dct(fCorrL);
FcorrR_new = dct(fCorrR);

VL = VL + (dt/2) * FcorrL_new;
VR = VR + (dt/2) * FcorrR_new;

W_solver = [UL; UR; VL; VR];

% --- Compare ---
rel_error = norm(W_matrix - W_solver) / norm(W_solver);
fprintf('Relative Error between G matrix and solver: %e\n', rel_error);

% Also verify dct(idct(x)) == x
x = randn(NxSub, 1);
fprintf('dct(idct(x)) error: %e\n', norm(dct(idct(x)) - x));
