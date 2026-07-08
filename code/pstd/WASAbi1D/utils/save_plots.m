function save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v)

    % Extracting test case data
    test_case = test_case_data.test_case;
    alpha_abs = test_case_data.alpha_abs;

    % Extracting simulation parameters
    merge_left = simulation_parameters.merge_left;
    merge_right = simulation_parameters.merge_right;
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    DD = simulation_parameters.DD;
    space_order = simulation_parameters.space_order;
    
    damped = alpha_abs ~= 0;

    % Save figures as images
    if DD == false
        merge_left_str = 'noDD';
        merge_right_str = 'noDD';
    else
        if merge_left == 1
            merge_left_str = 'pre_merge';
        else
            merge_left_str = 'post_merge';
        end

        if merge_right == 1
            merge_right_str = 'pre_merge';
        else
            merge_right_str = 'post_merge';
        end
    end
    
    if method_left == 1
        method_left_str = 'FDTD_2ord';
    elseif method_left == 2
        method_left_str = 'FDTD_1ord';
    elseif method_left == 3
        method_left_str = 'Fourier_2ord';
    elseif method_left == 4
        method_left_str = 'Fourier_1ord';
    elseif method_left == 5
        method_left_str = 'PML';
    end
    
    if method_right == 1
        method_right_str = 'FDTD_2ord';
    elseif method_right == 2
        method_right_str = 'FDTD_1ord';
    elseif method_right == 3
        method_right_str = 'Fourier_2ord';
    elseif method_right == 4
        method_right_str = 'Fourier_1ord';
    elseif method_right == 5
        method_right_str = 'PML';
    end
    
    if damped == true
        damping_str = 'damped';
    else
        damping_str = 'undamped';
    end
    
    % Create folder with current filename
    foldername = sprintf('test=%s__dh=%.4f_dt=%.4f__%s__left=%s_%s__right=%s_%s__space_order=%d', num2str(test_case), dh, dt, damping_str, method_left_str, merge_left_str, method_right_str, merge_right_str, space_order);

    % Create images folder if not exists
    if ~exist('images', 'dir')
        mkdir('images')
    end
    
    mkdir(fullfile('images', foldername));
    
    % Save figures as images in the folder
    filename_p_top = fullfile('images', foldername, 'pressure_top.png');
    filename_p_3d_a = fullfile('images', foldername, 'pressure_3d_a.png');
    filename_p_3d_b = fullfile('images', foldername, 'pressure_3d_b.png');
    filename_v_top = fullfile('images', foldername, 'velocity_top.png');
    filename_v_3d_a = fullfile('images', foldername, 'velocity_3d_a.png');
    filename_v_3d_b = fullfile('images', foldername, 'velocity_3d_b.png');
        
    view_top = true;
    view_3d = true;
    
    if view_top
        % Top view for pressure
        figure(fig_p)
        view(0,90)
        saveas(gcf, filename_p_top);
        
        % Top view for velocity
        figure(fig_v)
        view(0,90)
        saveas(gcf, filename_v_top);
    end
    
    if view_3d
        % 3D view A for pressure
        figure(fig_p)
        view(3)
        saveas(fig_p, filename_p_3d_a);
        
        % 3D view A for velocity
        figure(fig_v)
        view(3)
        saveas(gcf, filename_v_3d_a);


        % 3D view B for pressure
        figure(fig_p)
        view(30,30)
        saveas(fig_p, filename_p_3d_b);
        
        % 3D view B for velocity
        figure(fig_v)
        view(30,30)
        saveas(gcf, filename_v_3d_b);
    end

end