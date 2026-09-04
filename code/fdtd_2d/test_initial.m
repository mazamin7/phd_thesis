addpath('src');
Lx = pi; Ly = pi; Nx = 30; Ny = 30; dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
c = 1.0; gamma = 0.5; nu = 0.1; alpha = (gamma + 2*nu)/2; omega_mode = sqrt(2*c^2 - alpha^2);
dt = 0.5 * min(dx, dy) / c; T_end = 2.0;
Nt = floor(T_end / dt) + 1; dt = T_end / (Nt - 1);
u0_fun = @(X, Y) exp(-alpha*0) .* cos(omega_mode*0) .* sin(X) .* sin(Y);
v0_fun = @(X, Y) (-alpha*exp(-alpha*0) .* cos(omega_mode*0) - omega_mode*exp(-alpha*0) .* sin(omega_mode*0)) .* sin(X) .* sin(Y);
f_fun = @(X, Y, t) zeros(size(X));
bcType = 'dirichlet'; delta = 2 * dx;

[~, ~, t_grid, u_mono, ~] = solver_2d(u0_fun, v0_fun, f_fun, dx, dy, dt, Lx, Ly, T_end, c, gamma, nu, 0, 0, bcType, bcType, bcType, bcType, @(t)0, @(t)0, @(t)0, @(t)0);

mid = round((Lx/2)/dx) * dx; delta_half_L = ceil( (delta/2)/dx ) * dx; delta_half_R = floor( (delta/2)/dx ) * dx;
xL0 = 0; xL1 = mid + delta_half_L; xR0 = mid - delta_half_R; xR1 = Lx;
idxL_R = round((xL1 - xL0)/dx) + 1; idxR_L = round((xR0 - xL0)/dx) + 1; idxR_R = Nx;
uL = u_mono(:, 1:idxL_R, :); uR = u_mono(:, idxR_L:idxR_R, :);
idxL = round((xR0 - xL0)/dx) + 1; q = 0.1283; r = 1.5;

robinLeftVals = zeros(Ny, Nt);
for n = 3:Nt
    uxL_neg = (3*uL(:, idxL, n) - 4*uL(:, idxL+1, n) + uL(:, idxL+2, n)) / (2*dx);
    utL = (3*uL(:, idxL, n) - 4*uL(:, idxL, n-1) + uL(:, idxL, n-2)) / (2*dt);
    robinLeftVals(:, n) = uxL_neg + q*utL + r*uL(:, idxL, n);
end
robinLeftVals(:, 1:2) = repmat(robinLeftVals(:, 3), 1, 2);
gLR = @(t) interp1(t_grid, robinLeftVals', t, 'linear', 'extrap')';

[~, ~, ~, uR_new, ~] = solver_2d( ...
    @(X_sub, Y_sub) u0_fun(X_sub + xR0, Y_sub), ...
    @(X_sub, Y_sub) v0_fun(X_sub + xR0, Y_sub), ...
    @(X_sub, Y_sub, t) f_fun(X_sub + xR0, Y_sub, t), ...
    dx, dy, dt, xR1-xR0, Ly, T_end, ...
    c, gamma, nu, q, r, ...
    'robinbdf', 'dirichlet', 'dirichlet', 'dirichlet', ...
    gLR, @(t) zeros(Ny,1), @(t) zeros(idxR_R-idxR_L+1,1), @(t) zeros(idxR_R-idxR_L+1,1));

Nx_sub = idxR_R - idxR_L + 1;
[A, B, C] = interior_assembly_2d(Nx_sub, Ny, dx, dy, dt, c, gamma, nu);
[A_mod, B_mod, C_mod, rhsFunBc] = apply_bc_2d( ...
    A, B, C, Nx_sub, Ny, dx, dy, dt, q, r, ...
    'robinbdf', 'dirichlet', 'dirichlet', 'dirichlet', ...
    gLR, @(t) zeros(Ny,1), @(t) zeros(Nx_sub,1), @(t) zeros(Nx_sub,1));

uCurr = uR_new(:,:,2); uCurr = uCurr(:);
uPrev = uR_new(:,:,1); uPrev = uPrev(:);
rhs = B_mod * uCurr + C_mod * uPrev + rhsFunBc(t_grid(2));
u3_manual = A_mod \ rhs;
u3_manual = reshape(u3_manual, Ny, Nx_sub);

diff = abs(uR_new(:,:,3) - u3_manual);
[max_res, max_idx] = max(diff(:));
fprintf('Max difference between manual and solver_2d: %.4e\n', max_res);
