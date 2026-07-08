% compute_cfl.m
clear all; close all; clc;
addpath('WASAbi1D/utils');

orders = [2, 4, 6, 8];
cfl_limits = zeros(size(orders));

Nx = 100;
L = 1;
dx = L/Nx;
c = 1;

gamma = 0;
nu = 0;

NxSub = Nx/2;
Lsub = L/2;
km = (0:NxSub-1)' * pi / Lsub;
lambda = -(km.^2); % Since lambda is the Laplacian eigenvalue

% State vector ordering: W = [UL; UR; VL; VR]
% Size is 4*NxSub = 2*Nx
state_size = 2*Nx;

for k = 1:length(orders)
    space_order = orders(k);
    C_res = get_residue_matrix(Nx, space_order);
    
    % Construct the modal interface operator A (Nx by Nx)
    A = zeros(Nx, Nx);
    
    for col = 1:Nx % We loop over U part to find the mapping to Fcorr
        U_test = zeros(Nx, 1);
        U_test(col) = 1; % Impulse on U
        
        UL = U_test(1:NxSub);
        UR = U_test(NxSub+1:Nx);
        
        % Exact solver sequence
        uL = idct(UL);
        uR = idct(UR);
        
        u = [uL; uR];
        
        f = (c^2/dx^2) * (C_res * u);
        
        FL = dct(f(1:NxSub));
        FR = dct(f(NxSub+1:end));
        
        A(1:NxSub, col) = FL;
        A(NxSub+1:end, col) = FR;
    end
    
    % Bisection search for max CFL
    cfl_low = 0.0;
    cfl_high = 2.0; % Upper bound
    max_cfl = 0;
    
    for iter = 1:50
        CFL = (cfl_low + cfl_high) / 2;
        dt = CFL * dx / c;
        
        % Construct global drift matrix D directly from exact propagators
        [S_cells, ~] = modal_propagators(lambda, dt, c, gamma, nu);
        
        D = zeros(state_size, state_size);
        
        % The state is [UL; UR; VL; VR]
        % For mode j in left domain: UL(j) is index j, VL(j) is index Nx+j
        % For mode j in right domain: UR(j) is index NxSub+j, VR(j) is index Nx+NxSub+j
        
        for j = 1:NxSub
            S_mat = S_cells{j};
            
            % Left domain mapping:
            D(j, j) = S_mat(1,1);
            D(j, Nx+j) = S_mat(1,2);
            D(Nx+j, j) = S_mat(2,1);
            D(Nx+j, Nx+j) = S_mat(2,2);
            
            % Right domain mapping:
            D(NxSub+j, NxSub+j) = S_mat(1,1);
            D(NxSub+j, Nx+NxSub+j) = S_mat(1,2);
            D(Nx+NxSub+j, NxSub+j) = S_mat(2,1);
            D(Nx+NxSub+j, Nx+NxSub+j) = S_mat(2,2);
        end
        
        % Kick matrix constructed explicitly as block matrix
        K = [eye(Nx),   zeros(Nx);
             (dt/2)*A,  eye(Nx)];
        
        % Amplification matrix exactly matching ard_solver.m sequence
        G = K * D * K;
        
        % Use eig for stability
        rho = max(abs(eig(full(G))));
        
        % We use 1e-7 threshold because of floating point noise for multiple eigenvalues at 1.0
        if rho <= 1.0000001
            max_cfl = CFL;
            cfl_low = CFL;
        else
            cfl_high = CFL;
        end
    end
    
    cfl_limits(k) = max_cfl;
    fprintf('Space Order %d: Exact CFL Limit = %.4f\n', space_order, max_cfl);
end

disp(table(orders(:), cfl_limits(:), 'VariableNames', {'SpaceOrder', 'CFLLimit'}));