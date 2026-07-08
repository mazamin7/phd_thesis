function [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, diss, db_plot)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here

    addpath sim\FDTD\
    addpath sim\Fourier\

    % Extracting simulation parameters
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    space_order = simulation_parameters.space_order;

    % Extracting test case data
    len_x = test_case_data.len_x;
    len_t = test_case_data.len_t;
    c0 = test_case_data.c0;
    alpha1 = test_case_data.alpha1;
    alpha2 = test_case_data.alpha2;
    transmittivity = test_case_data.transmittivity;
    force_fun = test_case_data.force_fun;
    p_gt_fun = test_case_data.p_gt_fun;
    v_gt_fun = test_case_data.v_gt_fun;


    fourier_left = method_left == 2;
    fourier_right = method_right == 2;
    
    damped = alpha1 ~= 0;
    
        stable = check_stability(c0, dt, dh, space_order, fourier_left, diss);
        assert(stable, 'Stability condition not satisfied on the left');
    
        stable = check_stability(c0, dt, dh, space_order, fourier_right, diss);
        assert(stable, 'Stability condition not satisfied on the right');
    
    % Knowing simulation pars and test case, initialize simulation variables
    
    % Defining time and space axes
    N_t = floor(len_t / dt);
    N_x = floor(len_x / dh);
    
    x_axis = linspace(0+dh/2,len_x,N_x);
    t_axis = linspace(0,len_t,N_t);
    
    % Initialize solution, force, and boundary conditions data
    p = zeros(N_x,N_t+1);
    v = zeros(N_x,N_t+1);
    force = zeros(N_x,N_t+1);
    p_gt = zeros(N_x,N_t+1);
    v_gt = zeros(N_x,N_t+1);
    
    for n = 2:N_t+1
        force(:,n) = force_fun(x_axis, t_axis(n-1));
        p_gt(:,n) = p_gt_fun(x_axis, t_axis(n-1));
        v_gt(:,n) = v_gt_fun(x_axis, t_axis(n-1));
    end
    
    % Step 0
    % Imposing initial conditions (in position 1 we have stub)
    p(:,2) = p_gt(:,2);
    v(:,2) = v_gt(:,2);

        % Building residue matrix
        C = get_residue_matrix(N_x, space_order);
        
        % Initializing update methods
        if fourier_left == false
            data_left = init_FDTD(len_x/2, c0, dt, dh, space_order, alpha1, alpha2);
            disp("FDTD left")
        else
            data_left = init_Fourier(len_x/2, c0, dt, dh, alpha1, alpha2);
            disp("Fourier left")
        end
        
        if fourier_right == false
            data_right = init_FDTD(len_x/2, c0, dt, dh, space_order, alpha1, alpha2);
            disp("FDTD right")
        else
            data_right = init_Fourier(len_x/2, c0, dt, dh, alpha1, alpha2);
            disp("Fourier right")
        end

    if debug
        % Init figure
        f = figure();
        f.Position = [100, 100, 1200, 700];
    end


    % Steps 1:N_t-1 (shifted by 1 to simplify code)
    
    % Simulation loop
    for n = 2:N_t

            % Update force
            force_now = force(:,n);


            % Update pressure left
            if method_left == 1
                [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_FDTD(data_left, p(1:N_x/2,n), v(1:N_x/2,n), force_now(1:N_x/2));
            elseif method_left == 2
                [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_Fourier(data_left, p(1:N_x/2,n), v(1:N_x/2,n), force_now(1:N_x/2));
            end
            
            % Update pressure right
            if method_right == 1
                [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_FDTD(data_right, p(N_x/2+1:N_x,n), v(N_x/2+1:N_x,n), force_now(N_x/2+1:N_x));
            elseif method_right == 2
                [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_Fourier(data_right, p(N_x/2+1:N_x,n), v(N_x/2+1:N_x,n), force_now(N_x/2+1:N_x));
            end
            
            % Compute new residual
            residual = (c0 / dh)^2 * C * p(:,n+1);

            % Artificial dissipation for stability
            p(:,n+1) = diss * p(:,n+1);

            
            % Post-merge
            % using r^{n+1}, on v^{n+1}
            v(1:N_x/2,n+1) = v(1:N_x/2,n+1) + transmittivity^2 * dt * residual(1:N_x/2) / (1 + 2*dt*alpha1);
            v(N_x/2+1:N_x,n+1) = v(N_x/2+1:N_x,n+1) + transmittivity^2 * dt * residual(N_x/2+1:N_x) / (1 + 2*dt*alpha1);
        
        info_str = ['Instant [s]: ' num2str(n*dt, '%4.3f') ' / ' ...
                num2str(len_t, '%4.3f') ' ( ' num2str(n/N_t*100, '%4.1f') '% )'];

        if debug == true
            if mod(n-1,10) == 1
	            plot_snapshot(x_axis,len_x,p(:,n+1),v(:,n+1),c0,f,db_plot);
            end
            figure(f);
            sgtitle(info_str);
        else
            clc;
            disp(info_str);
        end
    
    end

    % Discarding stub time instant
    p = p(:,2:end);
    v = v(:,2:end);
    
end

