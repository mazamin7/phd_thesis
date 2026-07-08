function plot_ground_truth(test_case_data, db_plot)
    % Extracting test case data
    len_x = test_case_data.len_x;
    len_t = test_case_data.len_t;
    p_gt_fun = test_case_data.p_gt_fun;
    v_gt_fun = test_case_data.v_gt_fun;

    % Define the axes for plotting
    t_axis = linspace(0, len_t, 1e3);
    x_axis = linspace(0, len_x, 1e3);

    % Evaluate the functions on a grid
    [X, T] = meshgrid(x_axis, t_axis);
    p_gt = p_gt_fun(X, T)';
    v_gt = v_gt_fun(X, T)';

    plot_spacetime(t_axis, x_axis, p_gt, v_gt, db_plot, 'Ground Truth');
end
