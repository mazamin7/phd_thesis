addpath(genpath('../src'));
% Script to compute the optimal rho_opt over a 2D grid of (gamma, nu)

% Fixed parameters
c = 1;
dx = 0.05;
dt = 0.05;
L = 1;
T = 2;
delta = L/10;
a = (L-delta)/2;
b = a + delta;

% Grid definition
N_points = 30;
gamma_vals = linspace(0, 5, N_points);
nu_vals = logspace(-3, 1, N_points); % from 0.001 to 10

[Gamma, Nu] = meshgrid(gamma_vals, nu_vals);
RhoOpt_map = zeros(N_points, N_points);

fprintf('Starting 2D grid sweep over (gamma, nu)...\n');
for i = 1:N_points
    for j = 1:N_points
        [qOpt, rOpt, rhoOpt] = get_optimal_robin_params(c, Gamma(i,j), Nu(i,j), a, b, T, dx, dt);
        RhoOpt_map(i,j) = rhoOpt;
    end
    fprintf('Row %d/%d completed.\n', i, N_points);
end

fig = figure('Position', [100, 100, 600, 500]);
imagesc(gamma_vals, log10(nu_vals), RhoOpt_map);
set(gca, 'YDir', 'normal');
colormap jet;
colorbar;
xlabel('Telegrapher damping \gamma');
ylabel('log_{10}(\nu) [Viscoelastic damping]');
title('Optimal convergence factor \rho_{opt}');

% Format Y axis ticks nicely
yticks_log = -3:1;
yticks = yticks_log;
set(gca, 'YTick', yticks);

% Save figure
out_dir = '../../latex/Images/fdtd_snapshots/swr_parameters/';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
saveas(fig, fullfile(out_dir, 'swr_rho_opt_map.png'));
fprintf('Saved figure to %s\n', fullfile(out_dir, 'swr_rho_opt_map.png'));
