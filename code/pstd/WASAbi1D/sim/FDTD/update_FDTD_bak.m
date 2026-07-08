function [p_next, v_next] = update_FDTD(data, p_curr, v_curr, force_next)
    % UPDATE_FDTD Computes the next time step for pressure (p) and velocity (v) 
    % using a semi-implicit FDTD scheme, computing the Laplacian using an 
    % explicit matrix-vector product.
    
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
    center_idx = (kernel_len + 1) / 2; % e.g., 2 for a 3-point kernel, 4 for a 7-point kernel
    L_matrix = create_laplacian_matrix(laplacian_kernel(center_idx:end), N);
    
    if alpha2 ~= 0
        error("Viscoelastic damping FDTD not implemented. Matrix L for alpha2=0 is used.")
    end
    
    % --- 2. Pressure Update ---
    p_next = p_curr + dt * v_curr;
    
    % --- 3. Calculate Laplacian of p_next (Spatial term) ---
    % % Apply symmetric (Neumann) boundary conditions via padding
    % L_size = (length(laplacian_kernel) - 1) / 2;
    % p_next_symm = symmetrize(p_next, L_size);
    % 
    % % Calculate the discrete Laplacian via convolution
    % laplacian_p_next = conv(p_next_symm, laplacian_kernel, 'valid');

    laplacian_p_next = L_matrix * p_next;

    % figure(4)
    % hold off
    % plot(laplacian_p_next);
    % hold on
    % plot(laplacian_p_next_alt);

    % --- 4. Velocity Update ---
    v_next_numerator = v_curr + ...
                       (c^2 * dt / dh^2) * laplacian_p_next + ...
                       dt * force_next;
    
    telegrapher_denominator = 1 + 2 * dt * alpha1;
    v_next = v_next_numerator / telegrapher_denominator;
end

function p_symm = symmetrize(p, L_size)
% SYMMETRIZE Implements symmetric (Neumann) boundary conditions for convolution.
% Mirrors the data at the boundaries based on the convolution kernel size.
    
    % Mirror the boundary values
    left_boundary = p(1:L_size);
    right_boundary = p(end-L_size+1:end);
    
    % Concatenate mirrored boundaries and the core data
    p_symm = [flipud(left_boundary); p; flipud(right_boundary)];
end