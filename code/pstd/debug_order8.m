% debug_order8.m
clear all; close all; clc;
addpath('WASAbi1D/utils');

Nx = 100;
L = 1;
dx = L/Nx;
c = 1;
gamma = 0; nu = 0;
space_order = 8;
CFL = 0.9619;
dt = CFL * dx / c;

% 1. Compute G matrix
NxSub = Nx/2;
Lsub = L/2;
km = (0:NxSub-1)' * pi / Lsub;
lambda = -(km.^2);
state_size = 2*Nx;

C_res = get_residue_matrix(Nx, space_order);
A = zeros(Nx, Nx);
for col = 1:Nx
    U_test = zeros(Nx, 1);
    U_test(col) = 1;
    u = [idct(U_test(1:NxSub)); idct(U_test(NxSub+1:Nx))];
    f = (c^2/dx^2) * (C_res * u);
    A(1:NxSub, col) = dct(f(1:NxSub));
    A(NxSub+1:end, col) = dct(f(NxSub+1:end));
end

[S_cells, ~] = modal_propagators(lambda, dt, c, gamma, nu);
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
K = [eye(Nx), zeros(Nx); (dt/2)*A, eye(Nx)];
G = K * D * K;

eigs_all = eig(full(G));
rho = max(abs(eigs_all));
fprintf('Analytical rho for CFL=0.9619 is %.12f\n', rho);

% 2. Run ard_solver directly
T_sim = 25.0;
x0 = L/2;
sigma = 0.05;
u0_fun = @(x) exp(-((x-x0).^2)/(2*sigma^2));
v0_fun = @(x) zeros(size(x));
f_fun = @(x,t) zeros(size(x));

[~,~,~,~,u] = ard_solver(u0_fun, v0_fun, f_fun, dx, dt, L, T_sim, c, gamma, nu, 'neumann', 0, 0, space_order);
max_u = max(abs(u(:)));
fprintf('ard_solver Max Amplitude after T=25 is %e\n', max_u);
