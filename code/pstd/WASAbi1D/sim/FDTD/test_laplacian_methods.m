function test_laplacian_methods()
% TEST_LAPLACIAN_METHODS Compares the result of matrix multiplication (L*p) 
% with convolution (conv) for discrete Laplacian calculation with Neumann BCs.

    % Define common parameters
    N = 30; % Domain size
    dh = 1.0; % Spatial step size (for comparison scaling)
    
    % --- Input Vector: A test pressure field p_next ---
    % Must be a column vector of size N.
    p_test = (1:N)' * 0.1; 
    
    disp('======================================================');
    disp(['VERIFICATION TEST: N=', num2str(N), ', dh=', num2str(dh)]);
    disp('======================================================');

    % --- TEST 1: Standard 3-Point Laplacian Stencil ---
    disp('--- TEST 1: 3-Point Stencil (L_xx = [1 -2 1]/dh^2) ---');
    
    % 1a. Define Stencil (Full kernel for conv, Python format for matrix)
    KERNEL_3PT_FULL = [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90];
    LAPLACIAN_COEFFS_3PT = [-49/18, 3/2, -3/20, 1/90]; % Python format: [c_center, c_off1]

    % 1b. Matrix Method (L * p_test)
    L_matrix_3pt = create_laplacian_matrix(LAPLACIAN_COEFFS_3PT, N);
    
    % The matrix L contains the coefficients [1, -2, 1] without the 1/dh^2 scaling.
    % We must apply 1/dh^2 scaling here.
    matrix_result_3pt = (1 / dh^2) * (L_matrix_3pt * p_test);

    % 1c. Convolution Method (conv(symmetrize(p), kernel, 'valid'))
    L_size_3pt = (length(KERNEL_3PT_FULL) - 1) / 2; % L_size = 1
    p_symm_3pt = symmetrize_for_conv(p_test, L_size_3pt);
    
    % Apply 1/dh^2 scaling to the convolution result
    conv_result_3pt = (1 / dh^2) * conv(p_symm_3pt, KERNEL_3PT_FULL, 'valid');
    
    % 1d. Comparison
    max_diff_3pt = max(abs(matrix_result_3pt - conv_result_3pt));
    
    disp(['Matrix Result (Top 5):', newline, num2str(matrix_result_3pt(1:5)')]);
    disp(['Conv Result (Top 5):  ', newline, num2str(conv_result_3pt(1:5)')]);
    disp(['Max Difference (Should be near zero): ', num2str(max_diff_3pt)]);
    disp(' ');

    
    % --- TEST 2: High-Order 7-Point Laplacian Stencil ---
    disp('--- TEST 2: 7-Point Stencil (High-Order O(h^6)) ---');

    % 2a. Define Stencil (Full kernel for conv, Python format for matrix)
    c_center = -49 / 18; 
    c_1 = 3 / 2;         
    c_2 = -3 / 20;       
    c_3 = 1 / 90;        
    
    % Full 7-point kernel: [c3, c2, c1, c_center, c1, c2, c3]
    KERNEL_7PT_FULL = [c_3, c_2, c_1, c_center, c_1, c_2, c_3]; 
    LAPLACIAN_COEFFS_7PT = [c_center, c_1, c_2, c_3]; % Python format

    % 2b. Matrix Method (L * p_test)
    L_matrix_7pt = create_laplacian_matrix(LAPLACIAN_COEFFS_7PT, N);
    matrix_result_7pt = (1 / dh^2) * (L_matrix_7pt * p_test);

    % 2c. Convolution Method (conv(symmetrize(p), kernel, 'valid'))
    L_size_7pt = (length(KERNEL_7PT_FULL) - 1) / 2; % L_size = 3
    p_symm_7pt = symmetrize_for_conv(p_test, L_size_7pt);
    conv_result_7pt = (1 / dh^2) * conv(p_symm_7pt, KERNEL_7PT_FULL, 'valid');
    
    % 2d. Comparison
    max_diff_7pt = max(abs(matrix_result_7pt - conv_result_7pt));
    
    disp(['Matrix Result (Top 5):', newline, num2str(matrix_result_7pt(1:5)')]);
    disp(['Conv Result (Top 5):  ', newline, num2str(conv_result_7pt(1:5)')]);
    disp(['Max Difference (Should be near zero): ', num2str(max_diff_7pt)]);
    disp(' ');
end

% --- Auxiliary Function: Symmetric Padding for Convolution ---
function p_symm = symmetrize_for_conv(p, L_size)
% Implements symmetric (Neumann) boundary conditions for convolution.
    
    left_boundary = p(1:L_size);
    right_boundary = p(end-L_size+1:end);
    
    % Concatenate mirrored boundaries and the core data
    p_symm = [flipud(left_boundary); p; flipud(right_boundary)];
end

% --- Auxiliary Function: create_laplacian_matrix (Must be defined or accessible) ---
% NOTE: The correct function definition from the previous turn must be 
% present in your environment for this tester to run.
% The definition for create_laplacian_matrix is NOT included here for brevity.