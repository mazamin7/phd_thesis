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

[X, Y, t_grid, u_mono, v_mono] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, 0, 0, ...
    bcType, bcType, bcType, bcType, ...
    @(t)0, @(t)0, @(t)0, @(t)0);

mid = round((Lx/2)/dx) * dx;
delta_half_L = ceil( (delta/2)/dx ) * dx;
delta_half_R = floor( (delta/2)/dx ) * dx;
xL0 = 0; xL1 = mid + delta_half_L;
xR0 = mid - delta_half_R; xR1 = Lx;

idxL_L = 1; idxL_R = round((xL1 - xL0)/dx) + 1;
idxR_L = round((xR0 - xL0)/dx) + 1; idxR_R = Nx;

uL = u_mono(:, 1:idxL_R, :);
uR = u_mono(:, idxR_L:idxR_R, :);

idxL = round((xR0 - xL0)/dx) + 1;
idxR = round((xL1 - xR0)/dx) + 1;
Nt = length(t_grid);
q = 0.1283; r = 1.5;

robinLeftVals = zeros(Ny, Nt);
for n = 3:Nt
    uxL_neg = (3*uL(:, idxL, n) - 4*uL(:, idxL+1, n) + uL(:, idxL+2, n)) / (2*dx);
    utL = (3*uL(:, idxL, n) - 4*uL(:, idxL, n-1) + uL(:, idxL, n-2)) / (2*dt);
    robinLeftVals(:, n) = uxL_neg + q*utL + r*uL(:, idxL, n);
end
robinLeftVals(:, 1:2) = repmat(robinLeftVals(:, 3), 1, 2);

gLR = @(t) interp1(t_grid, robinLeftVals', t, 'linear', 'extrap')';

[~, ~, ~, uR_new, vR_new] = solver_2d( ...
    @(X_sub, Y_sub) u0_fun(X_sub + xR0, Y_sub), ...
    @(X_sub, Y_sub) v0_fun(X_sub + xR0, Y_sub), ...
    @(X_sub, Y_sub, t) f_fun(X_sub + xR0, Y_sub, t), ...
    dx, dy, dt, xR1-xR0, Ly, T_end, ...
    c, gamma, nu, q, r, ...
    'robinbdf', 'dirichlet', 'dirichlet', 'dirichlet', ...
    gLR, @(t) zeros(Ny,1), @(t) zeros(idxR_R-idxR_L+1,1), @(t) zeros(idxR_R-idxR_L+1,1));

err = max(max(abs(uR_new(:,:,end) - uR(:,:,end))));
fprintf('Fixed point error: %.4e\n', err);
