addpath('src');

Lx = pi; Ly = pi;
c = 1.0; gamma = 0.5; nu = 0.1;
alpha = (gamma + 2*nu)/2;
omega = sqrt(2*c^2 - alpha^2);

T_end = 2.0;

T_t = @(t) exp(-alpha*t) .* cos(omega*t);
T_prime_t = @(t) -alpha*exp(-alpha*t) .* cos(omega*t) - omega*exp(-alpha*t) .* sin(omega*t);

u0_fun = @(X, Y) T_t(0) .* sin(X) .* sin(Y);
v0_fun = @(X, Y) T_prime_t(0) .* sin(X) .* sin(Y);
f_fun = @(X, Y, t) zeros(size(X));
bcType = 'dirichlet';

N_values = [20, 40, 80, 160];
errors = zeros(size(N_values));
dh_values = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    Nx = N; Ny = N;
    dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
    dh_values(i) = dx;
    
    dt = 0.5 * dx / c;
    
    [X, Y, t_grid, u_array, v_array] = solver_2d( ...
        u0_fun, v0_fun, f_fun, ...
        dx, dy, dt, Lx, Ly, T_end, ...
        c, gamma, nu, ...
        bcType);

    u_exact = @(X, Y, t) T_t(t) .* sin(X) .* sin(Y);
    u_num = squeeze(u_array(:,:,end));
    u_ref = u_exact(X, Y, t_grid(end));
    
    errors(i) = max(abs(u_num(:) - u_ref(:)));
    fprintf("N=%d, dh=%.4f, dt=%.4f, Error=%e\n", N, dx, dt, errors(i));
end

rates = log2(errors(1:end-1) ./ errors(2:end));
for i = 1:length(rates)
    fprintf("Convergence rate from N=%d to N=%d: %.2f\n", N_values(i), N_values(i+1), rates(i));
end

figure;
loglog(dh_values, errors, '-o', 'LineWidth', 2);
hold on;
loglog(dh_values, dh_values.^2 * (errors(1)/dh_values(1)^2), '--k', 'LineWidth', 2);
xlabel('\Delta h');
ylabel('$L_\infty$ Error', 'Interpreter', 'latex');
legend('Numerical Error', '$\mathcal{O}(\Delta h^2)$ Reference', 'Interpreter', 'latex');
title('Convergence of 2D FDTD');
grid on;
saveas(gcf, 'convergence_rate.png');

