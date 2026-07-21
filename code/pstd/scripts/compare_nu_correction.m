% compare_nu_correction.m
clear all; close all; clc;
addpath(genpath('../src'));

%% Parameters
c = 1;
L = 1;
T = 2;
gamma = 0;
nu = 0.05;
bcType = 'neumannGhost';
q = 1/c;
r = 0;
space_order = 8; % Use 4th order interface

N = 50;
dx = L/N;

% Exact solution functions
uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
u0_fun = @(x) uExact(x,0);
v0_fun = @(x) zeros(size(x));
f_fun  = @(x,t) zeros(size(x));

% Time step
% The solver with explicit nu requires the parabolic limit
dt_parabolic = 0.1 * dx^2 / nu;
dt_hyperbolic = 0.9 * dx / c;
dt = min(dt_hyperbolic, dt_parabolic);

fprintf('Running with N=%d, dx=%.4f, dt=%.6f\n', N, dx, dt);

%% 1. Monolithic PSTD (Reference)
fprintf('Running monolithic PSTD...\n');
[x_grid, t_grid, u_mono, ~, U_mono, ~] = solver( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r);

%% 2. ARD With Exact Viscous Kick
fprintf('Running ARD with nu correction...\n');
[~, xL, xR, ~, u_ard_with, ~, ~, ~, ~] = ard_solver( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);

%% 3. ARD Without Viscous Kick (Omitted)
fprintf('Running ARD without nu correction...\n');
[~, ~, ~, ~, u_ard_without, ~, ~, ~, ~] = ard_solver_without_nu( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);

%% Compute Errors
u_mono_final = u_mono(:, end);
u_ard_with_final = u_ard_with(:, end);
u_ard_without_final = u_ard_without(:, end);

err_with = abs(u_ard_with_final - u_mono_final);
err_without = abs(u_ard_without_final - u_mono_final);

fprintf('\nMax Error (With Nu Kick): %.6e\n', max(err_with));
fprintf('Max Error (Without Nu Kick): %.6e\n', max(err_without));

%% Plot Comparison
h = figure('Position', [100 100 1000 400], 'Color', 'w');

subplot(1,2,1);
plot(x_grid, u_mono_final, 'k-', 'LineWidth', 2, 'DisplayName', 'Monolithic'); hold on;
plot(x_grid, u_ard_without_final, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ARD (No \nu Kick)');
plot(x_grid, u_ard_with_final, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'ARD (Exact \nu Kick)');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Wave Amplitude');
title(sprintf('Wavefield at T=%.1f', T));
legend('Location', 'best');

subplot(1,2,2);
semilogy(x_grid, err_without, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Error (No \nu Kick)'); hold on;
semilogy(x_grid, err_with, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Error (Exact \nu Kick)');
grid on;
xline(L/2, 'k:', 'Interface', 'HandleVisibility', 'off');
xlabel('x');
ylabel('Absolute Error |u_{ARD} - u_{mono}|');
title('Domain Decomposition Error');
legend('Location', 'best');

saveFolder = 'snapshots/nu_comparison';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
exportgraphics(h, fullfile(saveFolder, 'nu_kick_comparison.png'), 'Resolution', 300);

fprintf('Plot saved to %s\n', fullfile(saveFolder, 'nu_kick_comparison.png'));


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
