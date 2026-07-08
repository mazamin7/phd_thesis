% check_nx_trend.m
clear all; close all; clc;
addpath('WASAbi1D/utils');

orders = [2, 4, 6, 8];
Nxs = [20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 800];
L = 1;
c = 1;
gamma = 0; nu = 0;

results = zeros(length(orders), length(Nxs));

for n_idx = 1:length(Nxs)
    Nx = Nxs(n_idx);
    dx = L/Nx;
    NxSub = Nx/2;
    Lsub = L/2;
    km = (0:NxSub-1)' * pi / Lsub;
    lambda = -(km.^2);
    state_size = 2*Nx;
    
    for k = 1:length(orders)
        space_order = orders(k);
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
        
        cfl_low = 0.0;
        cfl_high = 2.0;
        max_cfl = 0;
        
        for iter = 1:30
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
            
            if rho <= 1.0000001
                max_cfl = CFL;
                cfl_low = CFL;
            else
                cfl_high = CFL;
            end
        end
        results(k, n_idx) = max_cfl;
        fprintf('Nx = %d, Order = %d, Limit = %.4f\n', Nx, space_order, max_cfl);
    end
end

var_names = cell(1, length(Nxs));
for n_idx = 1:length(Nxs)
    var_names{n_idx} = sprintf('Nx_%d', Nxs(n_idx));
end

disp('--- CFL Limits Asymptotic Trend ---');
disp(array2table(results, 'VariableNames', var_names, 'RowNames', {'Order_2', 'Order_4', 'Order_6', 'Order_8'}));

disp('--- Empirical Worst Case (Minimum) CFL Limit ---');
min_cfls = min(results, [], 2);
for k = 1:length(orders)
    fprintf('Order %d: %.4f\n', orders(k), min_cfls(k));
end
