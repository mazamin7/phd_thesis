addpath('src');

% Grid parameters
Lx = pi; Ly = pi;
Nx = 30; Ny = 30;
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
c = 1.0;
gamma = 0.5; nu = 1.0; % High viscoelastic damping

dt = 0.5 * min(dx, dy) / c;
T_end = 2.0;

% Multi-mode initial condition
u0_fun = @(X, Y) multi_mode(X, Y);
v0_fun = @(X, Y) zeros(size(X));
f_fun = @(X, Y, t) zeros(size(X));

% Subdomain setup
delta = 2 * dx;
bcType = 'dirichlet';

fprintf('Running monolithic...\n');
[~, ~, ~, u_mono, ~] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, 0, 0, ...
    bcType, bcType, bcType, bcType, ...
    @(t)0, @(t)0, @(t)0, @(t)0);

fprintf('Finding optimal params...\n');
qVals = linspace(0.1, 1.5/c, 100);
rVals = linspace(0, 1.5, 100);
Rho = zeros(length(rVals), length(qVals));

wmin = pi/T_end;
wmax = pi/dt;
omega_vals = logspace(log10(wmin), log10(wmax), 200);
ky_vals = (1:Ny) * pi / Ly;
[W, KY] = meshgrid(omega_vals, ky_vals);
k2d = sqrt((W.^2 - 1i*gamma*W)./(c^2 + 1i*nu*W) - KY.^2);

for i = 1:length(rVals)
    for j = 1:length(qVals)
        r = rVals(i); q = qVals(j);
        z = q * (1i*W) + r;
        rho_val = abs((z - 1i*k2d) ./ (z + 1i*k2d) .* exp(-1i*k2d*delta));
        Rho(i, j) = max(rho_val(:));
    end
end
[min_rho, idx] = min(Rho(:));
[r_idx, q_idx] = ind2sub(size(Rho), idx);
qOpt = qVals(q_idx); rOpt = rVals(r_idx);
fprintf('Optimum: q = %.4f, r = %.4f (rho = %.4f)\n', qOpt, rOpt, min_rho);

iterations = 15;
qClassic = 1/c; rClassic = 0;

fprintf('Running SWR (classical)...\n');
[~, ~, ~, u_iter_classic, ~] = swr_solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, ...
    qClassic, rClassic, ...
    delta, iterations);

fprintf('Running SWR (optimized)...\n');
[~, ~, ~, u_iter_opt, ~] = swr_solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, ...
    qOpt, rOpt, ...
    delta, iterations);

errClassic = zeros(iterations, 1);
errOpt = zeros(iterations, 1);

for k = 1:iterations
    diffClassic = u_iter_classic(:, :, end, k) - u_mono(:, :, end);
    errClassic(k) = max(abs(diffClassic(:)));
    
    diffOpt = u_iter_opt(:, :, end, k) - u_mono(:, :, end);
    errOpt(k) = max(abs(diffOpt(:)));
end

figure('Color', 'w');
semilogy(1:iterations, errClassic, 'o-', 'LineWidth', 2);
hold on;
semilogy(1:iterations, errOpt, 's-', 'LineWidth', 2);
grid on; box on;
xlabel('Schwarz iteration');
ylabel('L_\infty error');
legend('Robin q=1/c, r=0', 'Optimized Robin', 'Location', 'northeast');
title('2D SWR convergence');

saveas(gcf, 'swr_convergence_comparison_2d.png');
fprintf('Done.\n');

function u = multi_mode(X, Y)
    u = zeros(size(X));
    for m = 1:2:9
        for n = 1:2:9
            u = u + (1/(m^2*n^2)) * sin(m*X) .* sin(n*Y);
        end
    end
end
