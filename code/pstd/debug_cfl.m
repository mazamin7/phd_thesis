clear all; close all; clc;
addpath('WASAbi1D/utils');

Nx = 100;
L = 1;
dx = L/Nx;
c = 1;
gamma = 0;
nu = 0;
NxSub = Nx/2;
Lsub = L/2;
km = (0:NxSub-1)' * pi / Lsub;
lambda = -(km.^2); 
state_size = 2*Nx;

space_order = 8;
C_res = get_residue_matrix(Nx, space_order);
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

cfl_low = 0.0;
cfl_high = 2.0; 
for iter = 1:20
    CFL = (cfl_low + cfl_high) / 2;
    dt = CFL * dx / c;
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
    rho = max(abs(eig(full(G))));
    fprintf('iter %d: CFL=%.4f rho=%.10f\n', iter, CFL, rho);
    if rho <= 1.000000001
        cfl_low = CFL;
    else
        cfl_high = CFL;
    end
end
