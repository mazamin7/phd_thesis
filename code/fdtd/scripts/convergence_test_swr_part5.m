addpath(genpath('../src'));
% Script to compute the mapping from (gamma, nu) to optimal (q,r)

% Fixed parameters
c = 1;
dx = 0.05;
dt = 0.05;
% dt = 1e-3;
L = 1;
T = 2;
delta = L/10;
a = (L-delta)/2;
b = a + delta;

% Grid definition
N_gamma = 30;
N_nu = 30;
% gamma_vals = linspace(0, 5, N_gamma);
gamma_vals = logspace(0, 1, N_gamma);
nu_vals = logspace(-1, 1, N_nu);

[Gamma, Nu] = meshgrid(gamma_vals, nu_vals);
QOpt_map = zeros(N_nu, N_gamma);
ROpt_map = zeros(N_nu, N_gamma);

fprintf('Starting 2D grid sweep over (gamma, nu) for Q, R mapping...\n');
for i = 1:N_nu
    for j = 1:N_gamma
        if i == 1 && j == 1
            x0_guess = [];
        elseif i == 1
            x0_guess = [QOpt_map(i, j-1), ROpt_map(i, j-1)];
        else
            x0_guess = [QOpt_map(i-1, j), ROpt_map(i-1, j)];
        end
        
        [qOpt, rOpt, ~] = get_optimal_robin_params(c, Gamma(i,j), Nu(i,j), a, b, T, dx, dt, x0_guess);
        QOpt_map(i,j) = qOpt;
        ROpt_map(i,j) = rOpt;
    end
    fprintf('Row %d/%d completed.\n', i, N_nu);
end

out_dir = '../../../latex/Images/fdtd_snapshots/swr_parameters/';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Plot 1: Scatter of all (q, r) points
fig1 = figure('Position', [100, 100, 600, 500]);
scatter(QOpt_map(:), ROpt_map(:), 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Optimal q');
ylabel('Optimal r');
title('Mapping of all (\gamma, \nu) into (q, r) space');
grid on;
saveas(fig1, fullfile(out_dir, 'swr_qr_mapping_all.png'));
fprintf('Saved figure to %s\n', fullfile(out_dir, 'swr_qr_mapping_all.png'));

%% Plot 2: Coarse gamma, color mapped by nu
% Select 6 indices of gamma roughly evenly spaced
gamma_idx = round(linspace(1, N_gamma, 10));

fig2 = figure('Position', [200, 150, 700, 500]);
hold on;
cmap = jet(N_nu);

for g_i = gamma_idx
    % The curve for a fixed gamma over all nu
    q_curve = QOpt_map(:, g_i);
    r_curve = ROpt_map(:, g_i);
    
    % Draw the line to connect the dots
    plot(q_curve, r_curve, 'k-', 'LineWidth', 0.5);
    
    % Draw the scatter points with color based on nu
    scatter(q_curve, r_curve, 30, log10(nu_vals'), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % Add text label at the end of the curve (highest nu)
    text(q_curve(end), r_curve(end), sprintf('  \\gamma = %.1f', gamma_vals(g_i)), 'FontSize', 10);
end

colormap jet;
cb = colorbar;
ylabel(cb, 'log_{10}(\nu) [Viscoelastic damping]');
xlabel('Optimal q');
ylabel('Optimal r');
title('Optimal (q,r) trajectories for constant \gamma');
grid on;
saveas(fig2, fullfile(out_dir, 'swr_qr_mapping_gamma_curves.png'));
fprintf('Saved figure to %s\n', fullfile(out_dir, 'swr_qr_mapping_gamma_curves.png'));

%% Plot 3: Coarse nu, color mapped by gamma
% Select 5 indices of nu roughly evenly spaced in log scale
nu_idx = round(linspace(1, N_nu, 10));

fig3 = figure('Position', [300, 200, 700, 500]);
hold on;

for n_i = nu_idx
    % The curve for a fixed nu over all gamma
    q_curve = QOpt_map(n_i, :);
    r_curve = ROpt_map(n_i, :);
    
    % Draw the line to connect the dots
    plot(q_curve, r_curve, 'k-', 'LineWidth', 0.5);
    
    % Draw the scatter points with color based on gamma
    scatter(q_curve, r_curve, 30, gamma_vals, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % Add text label at the start of the curve (lowest gamma)
    text(q_curve(1), r_curve(1), sprintf(' \\nu = 10^{%.1f}', log10(nu_vals(n_i))), 'FontSize', 10);
end

colormap jet;
cb = colorbar;
ylabel(cb, 'Telegrapher damping \gamma');
xlabel('Optimal q');
ylabel('Optimal r');
title('Optimal (q,r) trajectories for constant \nu');
grid on;
saveas(fig3, fullfile(out_dir, 'swr_qr_mapping_nu_curves.png'));
fprintf('Saved figure to %s\n', fullfile(out_dir, 'swr_qr_mapping_nu_curves.png'));

fprintf('All mappings completed.\n');
