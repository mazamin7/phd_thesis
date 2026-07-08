function L = create_laplacian_matrix(python_coeffs, N)
% CREATE_LAPLACIAN_MATRIX Constructs the sparse Laplacian matrix L 
% with Neumann BCs enforced using the boundary reflection method.
% NOTE: Input python_coeffs MUST be in the format [c_center, c_off1, c_off2, ...]

    size_N = N;
    
    % --- CRITICAL FIX 1: Use the Python half-width definition and remove strict odd-length check. ---
    % python_coeffs length must be >= 1.
    L_c = fix(length(python_coeffs) - 1); 
    
    if L_c < 0 
        error('Invalid Stencil: Input coeffs must contain at least the center coefficient [c_center].');
    end
    
    row_size = size_N + 2 * L_c; 

    % --- 1. Create the base row for the circulant matrix ---
    row = zeros(1, row_size); 
    
    % Center coefficient is coeffs(1) in MATLAB / coeffs[0] in Python
    row(L_c + 1) = python_coeffs(1); 
    for i = 1:L_c
        % Off-diagonal coefficients are coeffs(i + 1) in MATLAB / coeffs[i] in Python
        row(L_c + 1 + i) = python_coeffs(i + 1); 
        row(L_c + 1 - i) = python_coeffs(i + 1); 
    end

    % --- 2. Build the full (N x N+2*L_c) matrix (mat_temp) ---
    mat_temp = zeros(size_N, row_size);
    for n = 1:size_N
        % Literal Python translation: torch.roll(row, n) -> circshift(row, [0, n-1])
        mat_temp(n, :) = circshift(row, [0, n - 1]); 
    end

    % --- 3. Extract the N x N central matrix (mat_N) ---
    mat_N = mat_temp(:, L_c + 1 : L_c + size_N);
    mat_B = zeros(size_N, size_N); % The boundary correction matrix

    % --- 4. Reflection Block Slicing and Addition (mat_B) ---
    if L_c > 0
        % Left Reflection (Top-Left Corner)
        left_overflow = mat_temp(1:L_c, 1:L_c);
        mat_B(1:L_c, 1:L_c) = mat_B(1:L_c, 1:L_c) + fliplr(left_overflow); 

        % Right Reflection (Bottom-Right Corner)
        right_overflow = mat_temp(size_N-L_c+1 : size_N, size_N+L_c+1 : row_size);
        mat_B(size_N-L_c+1 : size_N, size_N-L_c+1 : size_N) = ...
            mat_B(size_N-L_c+1 : size_N, size_N-L_c+1 : size_N) + fliplr(right_overflow);
    end

    % --- 5. Final Matrix ---
    L = sparse(mat_N + mat_B); 
end