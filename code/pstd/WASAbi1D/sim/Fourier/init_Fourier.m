function Fourier_data = init_Fourier(len_x, c0, dt, dh, air_alpha1, air_alpha2)
% Initializes the modal parameters for a 1D DCT-based acoustic solver.
% The parameters correspond to the 3D calculation in the C++ code, reduced to 1D.

    N = floor(len_x/dh); % Number of modes/grid points
    
    % The C++ code uses a 1-based index (k) for the modal wave number, 
    % where k=1 is the DC component (mode 0), and k=N is the highest mode.
    % In MATLAB, indices are 1 to N, corresponding to mode indices m=0 to N-1.
    
    % Modal index vector (m = 0, 1, ..., N-1)
    m = (0:N-1)'; 
    
    % Wave number terms in the C++ code for 1D: (k - 1) / lx2
    % where k-1 corresponds to m. lx2 is width^2 * dh^2.
    % The effective length L is len_x.
    lx2 = len_x^2;
    
    % 1D equivalent of the C++ w0 calculation:
    % w0 = c0 * pi * sqrt( (k-1)^2 / lx2_ )
    w0 = c0 * pi * m / len_x; 
    
    % Air absorption alpha
    % alpha = air_absorption_alpha1 + air_absorption_alpha2 * w0^2
    alpha = air_alpha1 + air_alpha2 * (w0.^2);
    
    % Pre-calculate constants for the time-stepping formulas
    inv_w02 = 1.0 ./ (w0.^2);
    
    % --- Initializing the main structure ---
    Fourier_data.N = N;
    Fourier_data.dt = dt;
    
    Fourier_data.w0 = w0;
    Fourier_data.alpha = alpha;
    Fourier_data.inv_w02 = inv_w02; 
    
    % Pre-calculate time stepping coefficients for each case
    
    % Case 1 & 2: idx == 0 (m = 0, DC component)
    % C++ separates by alpha == 0 or alpha > 0.
    % We will handle this separation inside the update function for clarity, 
    % but pre-calculate terms if alpha > 0.
    
    % Case 3: idx > 0 and alpha < w0 (Underdamped)
    w_under = sqrt(w0.^2 - alpha.^2);
    inv_w_under = 1.0 ./ w_under;
    alpha_sqr = alpha.^2;
    
    eatm_under = exp(-alpha .* dt);
    cwt_under = cos(w_under .* dt);
    swt_under = sin(w_under .* dt);
    
    Fourier_data.w_under = w_under;
    Fourier_data.inv_w_under = inv_w_under;
    Fourier_data.alpha_sqr = alpha_sqr;
    Fourier_data.eatm_under = eatm_under;
    Fourier_data.cwt_under = cwt_under;
    Fourier_data.swt_under = swt_under;
    
    % Case 4: idx > 0 and alpha > w0 (Overdamped)
    alphad = sqrt(alpha.^2 - w0.^2);
    alpha1 = alpha + alphad;
    alpha2 = alpha - alphad;
    
    eat1 = exp(-alpha1 .* dt);
    eat2 = exp(-alpha2 .* dt);
    
    Fourier_data.alphad = alphad;
    Fourier_data.alpha1 = alpha1;
    Fourier_data.alpha2 = alpha2;
    Fourier_data.eat1 = eat1;
    Fourier_data.eat2 = eat2;
end