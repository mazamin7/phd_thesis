% compare_strang_vs_exponential.m
clear all; close all; clc;
addpath(genpath('../src'));

%% Parameters (Experiment 4: Damped Trapezoid Wave)
c = 1;
L = 1;
T = 2;
gamma = 0;
nu = 0.05;
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

%% Run solver with Explicit Strang
fprintf('Running ARD with Explicit Strang Splitting...\n');
[x_grid, ~, ~, ~, u_exp, ~, ~, ~, ~] = ard_solver( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);
    
% Run solver with CN Strang
fprintf('Running ARD with CN Strang Splitting...\n');
[~, ~, ~, ~, u_cn, ~, ~, ~, ~] = ard_solver_cn_kick( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);

%% Compute Errors
u_mono_final = u_mono(:, end);
u_final_exp = u_exp(:, end);
u_final_cn = u_cn(:, end);

err_exp = max(abs(u_final_exp - u_mono_final));
err_cn = max(abs(u_final_cn - u_mono_final));

fprintf('\nMax Error (Explicit Strang Splitting): %e\n', err_exp);
fprintf('Max Error (CN Strang Splitting):       %e\n', err_cn);

%% Plot Comparison
h = figure('Position', [100 100 1000 400], 'Color', 'w');

subplot(1,2,1);
plot(x_grid, u_mono_final, 'k-', 'LineWidth', 2, 'DisplayName', 'Monolithic'); hold on;
plot(x_grid, u_final_exp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Explicit Strang');
plot(x_grid, u_final_cn, 'g:', 'LineWidth', 1.5, 'DisplayName', 'CN Strang');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Wave Amplitude');
title(sprintf('Wavefield at T=%.1f', T));
legend('Location', 'best');

subplot(1,2,2);
semilogy(x_grid, abs(u_final_exp - u_mono_final), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Error (Explicit Strang)'); hold on;
semilogy(x_grid, abs(u_final_cn - u_mono_final), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Error (CN Strang)');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Absolute Error |u_{ARD} - u_{mono}|');
title('Domain Decomposition Error');
legend('Location', 'best');

saveFolder = 'snapshots/strang_euler_vs_cn';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
exportgraphics(h, fullfile(saveFolder, 'strang_euler_vs_cn.png'), 'Resolution', 300);

fprintf('Plot saved to %s\n', fullfile(saveFolder, 'strang_euler_vs_cn.png'));


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
