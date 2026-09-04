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
iterations = 5;

[~, ~, ~, u_mono, ~] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, 0, 0, ...
    bcType, bcType, bcType, bcType, ...
    @(t)0, @(t)0, @(t)0, @(t)0);

[~, ~, ~, u_iter_opt, ~] = swr_solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, 0.1283, 1.5, ...
    delta, iterations);

for k = 1:iterations
    errOpt = max(max(abs(u_iter_opt(:,:,end,k) - u_mono(:,:,end))));
    fprintf('Iter %d: %.4e\n', k, errOpt);
end
