% compare_strang_vs_exponential.m
clear all; close all; clc;
addpath(genpath('../src'));

%% Parameters (Experiment 4: Damped Trapezoid Wave)
c = 1;
L = 1;
T = 2;
gamma = 0;
nu = 0;
bcType = 'neumannGhost';
q = 1/c;
r = 0;
space_order = 2;

N = 50;
dx = L/N;
% Use parabolic limit for stability!
dt_parabolic = 0.1 * dx^2 / nu;
dt_hyperbolic = 0.9 * dx / c;
dt = min(dt_hyperbolic, dt_parabolic);

fprintf('Running with N=%d, dx=%.4f, dt=%.6f\n', N, dx, dt);

% Exact solution functions
uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
u0_fun = @(x) uExact(x,0);
v0_fun = @(x) zeros(size(x));
f_fun  = @(x,t) zeros(size(x));

%% 1. Monolithic PSTD (Reference)
fprintf('Running monolithic PSTD...\n');
[x_grid, t_grid, u_mono, ~, U_mono, ~] = solver( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r);

%% 2. ARD With Strang Splitting
fprintf('Running ARD with Strang Splitting...\n');
[~, xL, xR, ~, u_strang, ~, ~, ~, ~] = ard_solver( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);

%% 3. ARD Without Splitting (Exponential Euler)
fprintf('Running ARD without splitting (Exponential Integrator)...\n');
[~, ~, ~, ~, u_no_split, ~, ~, ~, ~] = ard_solver_no_splitting( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);

%% Compute Errors
u_mono_final = u_mono(:, end);
u_strang_final = u_strang(:, end);
u_no_split_final = u_no_split(:, end);

err_strang = abs(u_strang_final - u_mono_final);
err_no_split = abs(u_no_split_final - u_mono_final);

fprintf('\nMax Error (Strang Splitting): %.6e\n', max(err_strang));
fprintf('Max Error (No Splitting / Exp. Euler): %.6e\n', max(err_no_split));

%% Plot Comparison
h = figure('Position', [100 100 1000 400], 'Color', 'w');

subplot(1,2,1);
plot(x_grid, u_mono_final, 'k-', 'LineWidth', 2, 'DisplayName', 'Monolithic'); hold on;
plot(x_grid, u_strang_final, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ARD (Strang Splitting)');
plot(x_grid, u_no_split_final, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'ARD (No Splitting)');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Wave Amplitude');
title(sprintf('Wavefield at T=%.1f', T));
legend('Location', 'best');

subplot(1,2,2);
semilogy(x_grid, err_strang, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Error (Strang Splitting)'); hold on;
semilogy(x_grid, err_no_split, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Error (No Splitting)');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Absolute Error |u_{ARD} - u_{mono}|');
title('Domain Decomposition Error');
legend('Location', 'best');

saveFolder = 'snapshots/splitting_comparison';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
exportgraphics(h, fullfile(saveFolder, 'splitting_comparison.png'), 'Resolution', 300);

fprintf('Plot saved to %s\n', fullfile(saveFolder, 'splitting_comparison.png'));


%% Helper Functions
function T_time = modal_time(k,c,gamma,nu)
    a = gamma + nu*k^2;
    disc = sqrt(a^2 - 4*c^2*k^2);
    lambda1 = (-a + disc)/2;
    lambda2 = (-a - disc)/2;
    A = -lambda2/(lambda1-lambda2);
    B =  lambda1/(lambda1-lambda2);
    T_time = @(t) A*exp(lambda1*t) + B*exp(lambda2*t);
end

function u = standing_trapezoid_neumann(x,t,M,L,c,gamma,nu)
    u = zeros(size(x));
    for j = 0:M
        m = 2*j + 1;
        k = m*pi/L;
        Tm = modal_time(k,c,gamma,nu);
        u = u + 16/pi^2 * cos(m*pi/4)/m^2 .* Tm(t) .* cos(k*x);
    end
end
