addpath('src');

Lx = pi; Ly = pi;
Nx = 50; Ny = 50;
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);

c = 1.0; gamma = 0.5; nu = 0.1;
alpha = (gamma + 2*nu)/2;
omega = sqrt(2*c^2 - alpha^2);

dt = 0.5 * min(dx, dy) / c;
T_end = 2.0;

T_t = @(t) exp(-alpha*t) .* cos(omega*t);
T_prime_t = @(t) -alpha*exp(-alpha*t) .* cos(omega*t) - omega*exp(-alpha*t) .* sin(omega*t);

u0_fun = @(X, Y) T_t(0) .* sin(X) .* sin(Y);
v0_fun = @(X, Y) T_prime_t(0) .* sin(X) .* sin(Y);
f_fun = @(X, Y, t) zeros(size(X));

bcType = 'dirichlet';

[X, Y, t_grid, u_array, v_array] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, ...
    bcType);

u_exact = @(X, Y, t) T_t(t) .* sin(X) .* sin(Y);

u_num = squeeze(u_array(:,:,end));
u_ref = u_exact(X, Y, t_grid(end));
err = max(abs(u_num(:) - u_ref(:)));
fprintf("Final refactored Max Error at T=%.2f: %e\n", t_grid(end), err);
