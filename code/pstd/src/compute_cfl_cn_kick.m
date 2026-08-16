% compute_cfl_exp_euler.m
clear all; close all; clc;
addpath(genpath('../src'));

orders = [2, 4, 6, 8];
cfl_hyperbolic = zeros(size(orders));
S_parabolic = zeros(size(orders));

Nx = 100;
L = 1;
dx = L/Nx;

NxSub = Nx/2;
Lsub = L/2;
km = (0:NxSub-1)' * pi / Lsub;
lambda = -(km.^2);
state_size = 2*Nx;

for k = 1:length(orders)
    space_order = orders(k);
    C_res = get_residue_matrix(Nx, space_order);
    
    A_base = zeros(Nx, Nx);
    for col = 1:Nx
        U_test = zeros(Nx, 1);
        U_test(col) = 1;
        
        uL = idct(U_test(1:NxSub));
        uR = idct(U_test(NxSub+1:Nx));
        u = [uL; uR];
        
        f = (1/dx^2) * (C_res * u);
        
        FL = dct(f(1:NxSub));
        FR = dct(f(NxSub+1:end));
        
        A_base(1:NxSub, col) = FL;
        A_base(NxSub+1:end, col) = FR;
    end
    
    %% Build global DCT/IDCT matrices
    % We need the IDCT operator Phi such that u_phys = Phi * U_modal
    % and DCT operator Phi_inv such that U_modal = Phi_inv * u_phys
    Phi = zeros(Nx, Nx);
    Phi_inv = zeros(Nx, Nx);
    for col = 1:NxSub
        vec = zeros(NxSub, 1); vec(col) = 1;
        Phi(1:NxSub, col) = idct(vec);
        Phi_inv(col, 1:NxSub) = dct(vec(1:NxSub));
        
        Phi(NxSub+1:end, NxSub+col) = idct(vec);
        Phi_inv(NxSub+col, NxSub+1:end) = dct(vec(1:NxSub));
    end
    
    %% HYPERBOLIC LIMIT
    c = 1; gamma = 0; nu = 0;
    
    cfl_low = 0.0; cfl_high = 2.0; max_cfl = 0;
    for iter = 1:50
        CFL = (cfl_low + cfl_high) / 2;
        dt = CFL * dx / c;
        [S_cells, T_cells] = modal_propagators(lambda, dt, c, gamma, nu);
        
        D = zeros(state_size, state_size);
        for j = 1:NxSub
            S_mat = S_cells{j};
            D(j, j) = S_mat(1,1); D(j, Nx+j) = S_mat(1,2);
            D(Nx+j, j) = S_mat(2,1); D(Nx+j, Nx+j) = S_mat(2,2);
            D(NxSub+j, NxSub+j) = S_mat(1,1); D(NxSub+j, Nx+NxSub+j) = S_mat(1,2);
            D(Nx+NxSub+j, NxSub+j) = S_mat(2,1); D(Nx+NxSub+j, Nx+NxSub+j) = S_mat(2,2);
        end
        
        dt_tau = dt/2;
        C_c2 = (c^2/dx^2) * C_res;
        C_nu = (nu/dx^2) * C_res;
        M_lhs = eye(Nx) - (dt_tau/2)*C_nu;
        M_rhs = eye(Nx) + (dt_tau/2)*C_nu;
        
        K21 = dt_tau * Phi_inv * (M_lhs \ C_c2) * Phi;
        K22 = Phi_inv * (M_lhs \ M_rhs) * Phi;
        
        K = [eye(Nx), zeros(Nx, Nx);
             K21,     K22];
        
        G = K * D * K;
        
        rho = max(abs(eig(full(G))));
        if rho <= 1.0000001
            max_cfl = CFL; cfl_low = CFL;
        else
            cfl_high = CFL;
        end
    end
    cfl_hyperbolic(k) = max_cfl;
    
    %% PARABOLIC LIMIT
    c = 1; gamma = 0; nu = 1;
    dx_test = 0.005;
    
    S_low = 0.0; S_high = 100.0; max_S = 0;
    for iter = 1:60
        S_num = (S_low + S_high) / 2;
        dt = S_num * dx_test^2 / nu;
        [S_cells, T_cells] = modal_propagators(lambda, dt, c, gamma, nu);
        
        D = zeros(state_size, state_size);
        for j = 1:NxSub
            S_mat = S_cells{j};
            D(j, j) = S_mat(1,1); D(j, Nx+j) = S_mat(1,2);
            D(Nx+j, j) = S_mat(2,1); D(Nx+j, Nx+j) = S_mat(2,2);
            D(NxSub+j, NxSub+j) = S_mat(1,1); D(NxSub+j, Nx+NxSub+j) = S_mat(1,2);
            D(Nx+NxSub+j, NxSub+j) = S_mat(2,1); D(Nx+NxSub+j, Nx+NxSub+j) = S_mat(2,2);
        end
        
        dt_tau = dt/2;
        C_c2 = (c^2/dx_test^2) * C_res;
        C_nu = (nu/dx_test^2) * C_res;
        M_lhs = eye(Nx) - (dt_tau/2)*C_nu;
        M_rhs = eye(Nx) + (dt_tau/2)*C_nu;
        
        K21 = dt_tau * Phi_inv * (M_lhs \ C_c2) * Phi;
        K22 = Phi_inv * (M_lhs \ M_rhs) * Phi;
        
        K = [eye(Nx), zeros(Nx, Nx);
             K21,     K22];
        
        G = K * D * K;
        
        rho = max(abs(eig(full(G))));
        if rho <= 1.0000001
            max_S = S_num; S_low = S_num;
        else
            S_high = S_num;
        end
    end
    S_parabolic(k) = max_S;
    
    fprintf('Space Order %d: CN-Strang CFL_max = %.4f, S_max = %.4f\n', space_order, max_cfl, max_S);
end
