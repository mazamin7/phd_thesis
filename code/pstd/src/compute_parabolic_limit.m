% compute_parabolic_limit.m
clear all; close all; clc;

orders = [2, 4, 6, 8];
S_limits = zeros(size(orders));

Nx = 200;
L = 1;
dx = L/Nx;
c = 1;
gamma = 0;
nu = 1; % Large viscosity to ensure parabolic limit strictly dominates

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
        
        uL = idct(U_test(1:NxSub));
        uR = idct(U_test(NxSub+1:Nx));
        u = [uL; uR];
        
        f = (c^2/dx^2) * (C_res * u);
        
        FL = dct(f(1:NxSub));
        FR = dct(f(NxSub+1:end));
        
        A(1:NxSub, col) = FL;
        A(NxSub+1:end, col) = FR;
    end
    
    % Bisection search for max S (parabolic dimensionless number)
    % S = nu * dt / dx^2
    S_low = 0.0;
    S_high = 10.0; % Increased upper bound
    max_S = 0;
    
    for iter = 1:50
        S = (S_low + S_high) / 2;
        dt = S * dx^2 / nu;
        
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
        
        B = (nu / c^2) * A;
        
        K = [eye(Nx),   zeros(Nx);
             (dt/2)*A,  eye(Nx) + (dt/2)*B];
             
        G = K * D * K;
        
        rho = max(abs(eig(full(G))));
        
        if rho <= 1.0000001
            max_S = S;
            S_low = S;
        else
            S_high = S;
        end
    end
    
    S_limits(k) = max_S;
    fprintf('Space Order %d: Exact Parabolic Limit S_max = %.4f\n', space_order, max_S);
end

disp(table(orders(:), S_limits(:), 'VariableNames', {'SpaceOrder', 'S_max'}));
