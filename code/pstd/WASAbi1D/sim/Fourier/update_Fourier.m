function [p_next, v_next] = update_Fourier(Fourier_data, p_curr, v_curr, force)
% Computes p_next and v_next using the modal update equations from C++.

    % Retrieving data
    N = Fourier_data.N;
    dt = Fourier_data.dt;
    w0 = Fourier_data.w0;
    alpha = Fourier_data.alpha;
    inv_w02 = Fourier_data.inv_w02;
    
    % Performing DCT (Discrete Cosine Transform type 2 used for the spatial part)
    DCT_type = 2;
    p_curr_modes = dct(p_curr,'Type',DCT_type);
    v_curr_modes = dct(v_curr,'Type',DCT_type);
    force_modes = dct(force,'Type',DCT_type); % Assuming force uses the same DCT type
    
    p_next_modes = zeros(N,1);
    v_next_modes = zeros(N,1);

    
    % --- Update solution in Modal Domain (Iterate through modes m = 0 to N-1) ---
    % MATLAB indices are 1 to N, corresponding to mode indices m=0 to N-1.

    % Mode m = 0 (idx = 0 in C++, MATLAB index 1)
    idx_dc = 1;
    
    if alpha(idx_dc) == 0.0
        % Case 1: idx == 0 and alpha == 0.0 (Undamped DC)
        % next_velocity_modes[0] = velocity_.modes_[0] + dt_ * force_.modes_[0];
        % next_pressure_modes[0] = dt_ * velocity_.modes_[0] + pressure_.modes_[0] + dt_*dt_/2 * force_.modes_[0];
        v_next_modes(idx_dc) = v_curr_modes(idx_dc) + dt * force_modes(idx_dc);
        p_next_modes(idx_dc) = dt * v_curr_modes(idx_dc) + p_curr_modes(idx_dc) + 0.5 * dt^2 * force_modes(idx_dc);
        
    elseif alpha(idx_dc) > 0.0
        % Case 2: idx == 0 and alpha > 0.0 (Damped DC)
        % e2at = exp(-2 * alpha * dt_);
        % next_velocity_modes[0] = e2at * velocity_.modes_[0] + (1 - e2at)/(2 * alpha) * force_.modes_[0];
        % next_pressure_modes[0] = (1 - e2at)/(2 * alpha) * velocity_.modes_[0] + pressure_.modes_[0] + ((e2at - 1) / (4 * alpha * alpha) + 1 / (2 * alpha) * dt_) * force_.modes_[0];
        
        alpha_dc = alpha(idx_dc);
        e2at = exp(-2 * alpha_dc * dt);
        inv_2a = 1.0 / (2 * alpha_dc);
        inv_4a2 = 1.0 / (4 * alpha_dc^2);
        
        v_next_modes(idx_dc) = e2at * v_curr_modes(idx_dc) + (1 - e2at) * inv_2a * force_modes(idx_dc);
        p_next_modes(idx_dc) = (1 - e2at) * inv_2a * v_curr_modes(idx_dc) + p_curr_modes(idx_dc) + ((e2at - 1) * inv_4a2 + inv_2a * dt) * force_modes(idx_dc);
        
    end

    % Modes m = 1 to N-1 (idx > 0 in C++, MATLAB indices 2 to N)
    idx_modes = 2:N; 
    
    % Pre-calculate separation criteria
    w0_modes = w0(idx_modes);
    alpha_modes = alpha(idx_modes);
    inv_w02_modes = inv_w02(idx_modes);
    
    % --- Case 3: idx > 0 and alpha < w0 (Underdamped/Oscillatory) ---
    underdamped_mask = alpha_modes < w0_modes;
    idx_under = idx_modes(underdamped_mask);
    
    if ~isempty(idx_under)
        % Retrieve pre-calculated coefficients
        inv_w_under = Fourier_data.inv_w_under(idx_under);
        alpha_sqr = Fourier_data.alpha_sqr(idx_under);
        eatm = Fourier_data.eatm_under(idx_under);
        cwt = Fourier_data.cwt_under(idx_under);
        swt = Fourier_data.swt_under(idx_under);
        w_modes = Fourier_data.w_under(idx_under);
        
        alpha_m = alpha_modes(underdamped_mask);
        inv_w02_m = inv_w02_modes(underdamped_mask);

        xe = force_modes(idx_under) .* inv_w02_m;
        
        % next_velocity_modes[idx] = eatm * (v_curr * (cwt - alpha * inv_w * swt) - (w + alpha_sqr * inv_w) * (p_curr - xe) * swt);
        v_next_modes(idx_under) = eatm .* (v_curr_modes(idx_under) .* (cwt - alpha_m .* inv_w_under .* swt) - ...
                                          (w_modes + alpha_sqr .* inv_w_under) .* (p_curr_modes(idx_under) - xe) .* swt);
        
        % next_pressure_modes[idx] = xe + eatm * ((p_curr - xe) * (cwt + alpha * inv_w * swt) + swt * inv_w * v_curr);
        p_next_modes(idx_under) = xe + eatm .* ((p_curr_modes(idx_under) - xe) .* (cwt + alpha_m .* inv_w_under .* swt) + ...
                                                swt .* inv_w_under .* v_curr_modes(idx_under));
    end

    % --- Case 4: idx > 0 and alpha > w0 (Overdamped/Non-oscillatory) ---
    overdamped_mask = alpha_modes > w0_modes;
    idx_over = idx_modes(overdamped_mask);
    
    if ~isempty(idx_over)
        % Retrieve pre-calculated coefficients
        alphad = Fourier_data.alphad(idx_over);
        alpha1 = Fourier_data.alpha1(idx_over);
        alpha2 = Fourier_data.alpha2(idx_over);
        eat1 = Fourier_data.eat1(idx_over);
        eat2 = Fourier_data.eat2(idx_over);
        
        inv_alphad = 1.0 ./ alphad;
        inv_w02_m = inv_w02_modes(overdamped_mask);

        % Note: C++ uses 0.5/alphad and 0.5 * ... to factor out the inv_w02 
        % in the velocity equation, so be careful with the scaling.
        
        % Terms used in both updates
        half_inv_alphad = 0.5 * inv_alphad;
        
        v_term1 = 0.5 * (eat1 + eat2);
        v_term2 = half_inv_alphad .* alpha(idx_over) .* (eat1 - eat2);
        v_term3 = -0.5 * (alpha1 .* eat1 + alpha2 .* eat2);
        v_term4 = -half_inv_alphad .* alpha(idx_over) .* (alpha2 .* eat2 - alpha1 .* eat1);
        v_term5 = inv_w02_m .* 0.5 .* (alpha1 .* eat1 + alpha2 .* eat2 + alpha(idx_over) .* inv_alphad .* (alpha2 .* eat2 - alpha1 .* eat1));

        % v_next_modes[idx] = (v_term1 + v_term2) * v_curr + (v_term3 + v_term4) * p_curr + v_term5 * force
        v_next_modes(idx_over) = (v_term1 + v_term2) .* v_curr_modes(idx_over) + ...
                                 (v_term3 + v_term4) .* p_curr_modes(idx_over) + ...
                                 v_term5 .* force_modes(idx_over);

        % p_next_modes[idx] = half_inv_alphad * (eat2 - eat1) * v_curr + 
        %                   (eat1 + eat2 + half_inv_alphad * (alpha2 * eat2 - alpha1 * eat1)) * p_curr + 
        %                   inv_w02 * (1 - alpha1 - alpha2 - half_inv_alphad * (alpha2 * eat2 - alpha1 * eat1)) * force 
        
        p_term1 = half_inv_alphad .* (eat2 - eat1);
        p_term2 = eat1 + eat2 + half_inv_alphad .* (alpha2 .* eat2 - alpha1 .* eat1);
        p_term3 = inv_w02_m .* (1 - alpha1 - alpha2 - half_inv_alphad .* (alpha2 .* eat2 - alpha1 .* eat1)); 
        
        p_next_modes(idx_over) = p_term1 .* v_curr_modes(idx_over) + ...
                                 p_term2 .* p_curr_modes(idx_over) + ...
                                 p_term3 .* force_modes(idx_over);
    end
    
    % Note: The C++ code doesn't explicitly handle the critical damping case (alpha == w0),
    % which would require a separate formula. For stability, the existing cases may
    % still be used with $\alpha$ slightly perturbed, but for a complete implementation, 
    % the critical damping case should also be implemented.
    % The given MATLAB only implements the underdamped case, we've fixed that.
    
    % --- Perform IDCT ---
    v_next = idct(v_next_modes,'Type',DCT_type);
    p_next = idct(p_next_modes,'Type',DCT_type);
end