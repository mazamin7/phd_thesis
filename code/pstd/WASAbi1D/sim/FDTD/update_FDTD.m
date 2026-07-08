function [p_next, v_next] = update_FDTD(data, p_curr, v_curr, force_next)
% UPDATE_FDTD Computes the next time step for pressure (p) and velocity (v) 
% using a STABLE IMPLICIT FDTD scheme, supporting both Telegrapher (alpha1) 
% and Viscoelastic (alpha2) damping via a matrix-based solution.

    % --- 1. Retrieve Simulation Parameters and Matrix ---
    laplacian_kernel = data.laplacian;
    c = data.c;
    dt = data.dt;
    dh = data.dh;
    alpha1 = data.alpha1;
    alpha2 = data.alpha2;
    N = length(p_curr);
    
    % --- Retrieve or Create Laplacian Matrix (L) ---
    kernel_len = length(laplacian_kernel);
    center_idx = (kernel_len + 1) / 2; 
    
    % L_matrix holds the unscaled spatial operator (L_xx) with Neumann BCs.
    % This call is assumed to be fixed and correct.
    L_matrix = create_laplacian_matrix(laplacian_kernel(center_idx:end), N);
    
    % --- 2. Pressure Update (Explicit step) ---
    % p_next = p_curr + dt * v_curr
    p_next = p_curr + dt * v_curr;
    
    % --- 3. Right-Hand Side (b) Construction (Explicit terms) ---
    
    % 3a. Calculate Laplacian of p_next (Spatial term)
    % laplacian_p_next = L_matrix * p_next (Unscaled coefficients applied to p_next)
    laplacian_p_next = L_matrix * p_next;
    
    % 3b. Build the RHS vector 'b'
    % b = v_curr + dt * [ (c^2 / dh^2) * L_xx(p_next) + f_next ]
    b = v_curr + dt * ( (c^2 / dh^2) * laplacian_p_next + force_next );

    
    % --- 4. Left-Hand Side (A) Construction (Implicit terms) ---
    
    % A corresponds to the implicit part of the equation: A * v_next = b
    % The terms contributing to A are:
    % 1. Identity (from v_next itself): I
    % 2. Telegrapher damping: 2 * dt * alpha1 * I
    % 3. Viscoelastic damping: - 2 * dt * c^2 * alpha2 / dh^2 * L_xx
    
    % Term 1 & 2: Main diagonal scalar component
    D_scalar = 1 + 2 * dt * alpha1;

    % Term 3: Viscoelastic Damping Operator Coefficient (scaled)
    % C_v_operator = - 2 * dt * c^2 * alpha2 / dh^2
    C_v_operator = - 2 * dt * c^2 * alpha2 / dh^2; 

    % The full matrix A is the sum of two operators:
    % A = D_scalar * I + C_v_operator * L_matrix
    
    % Note on Matrix Construction: The L_matrix is sparse, but D_scalar * I 
    % must be added to its diagonal.
    
    A = speye(N) * D_scalar + L_matrix * C_v_operator; 
    % A = (1 + 2*dt*alpha1)*I - (2*dt*c^2*alpha2/dh^2) * L_matrix
    
    % --- 5. Solve for v_next ---
    % The system is solved using the stable backslash operator.
    % If alpha2 = 0, L_matrix * C_v_operator becomes zero, and A = D_scalar * I.
    % The solution simplifies to v_next = b / D_scalar, which recovers the 
    % original explicit Telegrapher scheme (correctly handling the alpha2=0 case).
    
    % Solve A * v_next = b
    v_next = A \ b; 
    
    % --- Optional: Restore original non-implicit behavior if alpha2 is forced to zero ---
    % Since the implicit solution is stable and accurate even if alpha2=0, 
    % the 'if alpha2 ~= 0' check and error block are no longer needed.
end

function p_symm = symmetrize(p, L_size)
% SYMMETRIZE Implements symmetric (Neumann) boundary conditions for convolution.
    
    % Mirror the boundary values
    left_boundary = p(1:L_size);
    right_boundary = p(end-L_size+1:end);
    
    % Concatenate mirrored boundaries and the core data
    p_symm = [flipud(left_boundary); p; flipud(right_boundary)];
end