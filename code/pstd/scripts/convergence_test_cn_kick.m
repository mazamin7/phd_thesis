% convergence_test_exp_euler.m
clear all; close all; clc;
addpath(genpath('../src'));

%% Experiment Parameters
c = 1;
L = 1;
T = 0.5; % Short time to keep things fast
gamma = 0;
nu = 0.05;
bcType = 'neumannGhost';
q = 1/c;
r = 0;

N_array = [50, 100, 200, 400, 800];
dx_array = L ./ N_array;
% orders = [2, 4, 6, 8];
orders = [2];

errors = zeros(length(N_array), length(orders));

uExact = @(x,t) standing_trapezoid_neumann(x,t,500,L,c,gamma,nu);

for j = 1:length(orders)
    space_order = orders(j);
    fprintf('\n--- Space Order %d ---\n', space_order);
    
    for i = 1:length(N_array)
        N = N_array(i);
        dx = L/N;
        
        % Exponential Euler requires dt to satisfy parabolic limit
        dt = 0.5 * dx^2 / nu; % S = 0.5
        
        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun  = @(x,t) zeros(size(x));
        
        % Run CN Strang Splitting
        [x_grid, ~, ~, ~, u_ard, ~, ~, ~, ~] = ard_solver_cn_kick( ...
            u0_fun, v0_fun, f_fun, ...
            dx, dt, L, T, c, gamma, nu, bcType, q, r, space_order);
            
        u_final = u_ard(:, end);
        u_ref = uExact(x_grid, T);
        
        err = max(abs(u_final - u_ref));
        errors(i, j) = err;
        
        fprintf('N=%3d, dx=%.4f, dt=%.6f: Max Error = %.6e\n', N, dx, dt, err);
    end
end

%%
% Compute empirical convergence rates
rates = zeros(length(N_array)-1, length(orders));
for j = 1:length(orders)
    for i = 1:length(N_array)-1
        rates(i, j) = log2(errors(i,j) / errors(i+1,j));
    end
end

fprintf('\nEmpirical Convergence Rates:\n');
% disp(array2table(rates, 'VariableNames', {'Order2', 'Order4', 'Order6', 'Order8'}, ...
%     'RowNames', {'N:50->100', 'N:100->200', 'N:200->400', 'N:400->800'}));

% Plot
h = figure('Position', [100 100 600 500], 'Color', 'w');
loglog(dx_array, errors(:,1), 'o-', 'LineWidth', 2, 'DisplayName', 'Order 2'); hold on;
% loglog(dx_array, errors(:,2), 's-', 'LineWidth', 2, 'DisplayName', 'Order 4');
% loglog(dx_array, errors(:,3), '^-', 'LineWidth', 2, 'DisplayName', 'Order 6');
% loglog(dx_array, errors(:,4), 'd-', 'LineWidth', 2, 'DisplayName', 'Order 8');

% Reference lines
loglog(dx_array, 1e-1*dx_array.^2, 'k--', 'LineWidth', 1.5, 'DisplayName', '\mathcal{O}(\Delta x^2)');
loglog(dx_array, 1e-1*dx_array.^4, 'k-.', 'LineWidth', 1.5, 'DisplayName', '\mathcal{O}(\Delta x^4)');

grid on;
xlabel('\Delta x');
ylabel('L_\infty Error at T=0.5');
title('Spatial Convergence of CN Strang Splitting (Exp 4)');
legend('Location', 'southeast');

saveFolder = 'snapshots/convergence_cn_kick';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
exportgraphics(h, fullfile(saveFolder, 'convergence_rates.png'), 'Resolution', 300);
fprintf('Saved plot to %s\n', fullfile(saveFolder, 'convergence_rates.png'));

%% Helper
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
