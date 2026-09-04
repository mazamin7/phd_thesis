addpath('src');

Lx = pi; Ly = pi;
Nx = 30; Ny = 30;
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);

c = 1.0; gamma = 0.5; nu = 0.1;
alpha = (gamma + 2*nu)/2;
omega_mode = sqrt(2*c^2 - alpha^2);

dt = 0.5 * min(dx, dy) / c;
T_end = 2.0;

T_t = @(t) exp(-alpha*t) .* cos(omega_mode*t);
T_prime_t = @(t) -alpha*exp(-alpha*t) .* cos(omega_mode*t) - omega_mode*exp(-alpha*t) .* sin(omega_mode*t);

u0_fun = @(X, Y) T_t(0) .* sin(X) .* sin(Y);
v0_fun = @(X, Y) T_prime_t(0) .* sin(X) .* sin(Y);
f_fun = @(X, Y, t) zeros(size(X));

bcType = 'dirichlet';

delta = 2 * dx;
iterationsToMeasure = 5;

mid = Lx/2;
delta_half_L = ceil( (delta/2)/dx ) * dx;
delta_half_R = floor( (delta/2)/dx ) * dx;
a = mid - delta_half_R;
b = mid + delta_half_L;

% Monolithic reference
[X, Y, t_grid, u_mono, ~] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, 0, 0, ...
    bcType, bcType, bcType, bcType, ...
    @(t)0, @(t)0, @(t)0, @(t)0);

qVals = linspace(0, 1.5/c, 15);
rVals = linspace(0, 2, 15);

Err = zeros(length(rVals), length(qVals));
Rho = zeros(length(rVals), length(qVals));

wmin = pi/T_end;
wmax = pi/dt;
omega_vals = logspace(log10(wmin), log10(wmax), 200);
ky_vals = (1:Ny) * pi / Ly;

[W, KY] = meshgrid(omega_vals, ky_vals);
k2d = sqrt((W.^2 - 1i*gamma*W)./(c^2 + 1i*nu*W) - KY.^2);

for iq = 1:length(qVals)
    fprintf('%d / %d\n', iq, length(qVals));
    for ir = 1:length(rVals)
        q = qVals(iq);
        r = rVals(ir);

        Lambda = r + 1i*q*W;
        rho = abs((1i*k2d.*sinh(1i*k2d*a)-Lambda.*cosh(1i*k2d*a))./(1i*k2d.*sinh(1i*k2d*b)+Lambda.*cosh(1i*k2d*b))).^2 ...
             .*exp(-2*imag(k2d)*delta);
        Rho(ir, iq) = max(rho(:));

        [~, ~, ~, u_iter, ~] = swr_solver_2d( ...
            u0_fun, v0_fun, f_fun, ...
            dx, dy, dt, Lx, Ly, T_end, ...
            c, gamma, nu, q, r, ...
            delta, iterationsToMeasure);

        Err(ir, iq) = max(max(abs(u_iter(:,:,end,end) - u_mono(:,:,end))));
    end
end

[minErr,idx] = min(Err(:));
[irErr,iqErr] = ind2sub(size(Err),idx);
[minRho,idx] = min(Rho(:));
[irRho,iqRho] = ind2sub(size(Rho),idx);

fig1 = figure('Color','w');
imagesc(qVals, rVals, log10(Err));
axis xy; colorbar;
xlabel('q'); ylabel('r');
title(sprintf('log_{10} error after %d Schwarz iterations (2D)', iterationsToMeasure));
hold on;
plot(qVals(iqErr), rVals(irErr), 'wo', 'MarkerSize', 12, 'LineWidth', 2);
plot(qVals(iqRho), rVals(irRho), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
legend('Measured optimum', 'Continuous optimum');
exportgraphics(fig1, 'swr_empirical_error_2d.png', 'Resolution', 300);

fig2 = figure('Color','w');
imagesc(qVals, rVals, Rho);
axis xy; colorbar;
xlabel('q'); ylabel('r');
title('Predicted max convergence factor (2D)');
hold on;
plot(qVals(iqErr), rVals(irErr), 'wo', 'MarkerSize', 12, 'LineWidth', 2);
plot(qVals(iqRho), rVals(irRho), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
legend('Measured optimum', 'Continuous optimum');
exportgraphics(fig2, 'swr_predicted_rho_2d.png', 'Resolution', 300);
