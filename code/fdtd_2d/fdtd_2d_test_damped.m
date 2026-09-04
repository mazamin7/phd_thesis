Lx = pi; Ly = pi;
Nx = 50; Ny = 50;
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);

c = 1.0;
gamma = 0.5;
nu = 0.1;
alpha = (gamma + 2*nu)/2;
omega = sqrt(2*c^2 - alpha^2);

dt = 0.5 * min(dx, dy) / c;
T_end = 2.0;
Nt = floor(T_end / dt) + 1;
dt = T_end / (Nt - 1);

x = linspace(0, Lx, Nx)';
y = linspace(0, Ly, Ny)';
[X, Y] = meshgrid(x, y);

% Analytical solution
T_t = @(t) exp(-alpha*t) .* cos(omega*t);
T_prime_t = @(t) -alpha*exp(-alpha*t) .* cos(omega*t) - omega*exp(-alpha*t) .* sin(omega*t);
T_prime_prime_t = @(t) (alpha^2 - omega^2)*exp(-alpha*t) .* cos(omega*t) + 2*alpha*omega*exp(-alpha*t) .* sin(omega*t);

u_exact = @(X, Y, t) T_t(t) .* sin(X) .* sin(Y);
v_exact = @(X, Y, t) T_prime_t(t) .* sin(X) .* sin(Y);
utt_exact = @(X, Y, t) T_prime_prime_t(t) .* sin(X) .* sin(Y);

u0 = u_exact(X, Y, 0);
v0 = v_exact(X, Y, 0);
utt0 = utt_exact(X, Y, 0);

uCurr = u0(:);
uNext = uCurr + dt*v0(:) + 0.5*dt^2*utt0(:);
uPrev = uCurr;
uCurr = uNext;

% Matrices
Ix = speye(Nx); Iy = speye(Ny);
e_x = ones(Nx, 1); e_y = ones(Ny, 1);
Dxx = spdiags([e_x -2*e_x e_x], [-1 0 1], Nx, Nx) / dx^2;
Dyy = spdiags([e_y -2*e_y e_y], [-1 0 1], Ny, Ny) / dy^2;

% Note kron order: kron(A, B) does A \otimes B.
% If u is vectorized as u(:) where columns of U are stacked,
% Dxx applies to the rows of U (y-direction) or columns (x-direction)?
% meshgrid(x,y) produces X where rows vary along y, cols vary along x.
% So U(j, i) corresponds to y_j, x_i.
% When vectorizing u(:), it stacks columns. So the fast index is y, the slow index is x.
% Dyy acts on columns of U (fast index). Dxx acts on rows of U (slow index).
% Therefore D_total = kron(Ix, Dyy) + kron(Dxx, Iy). Wait.
% If fast index is y, then a 1D operator on y is kron(Ix, Dyy).
% A 1D operator on x is kron(Dxx, Iy).
% Let's verify: D_total = kron(Dxx, Iy) + kron(Ix, Dyy);
D = kron(Dxx, Iy) + kron(Ix, Dyy);
I = speye(Nx * Ny);

A = (1 + gamma*dt/2)*I - (nu*dt/2)*D;
B = 2*I + c^2*dt^2*D;
C = -(1 - gamma*dt/2)*I - (nu*dt/2)*D;

is_boundary = false(Ny, Nx);
is_boundary(1, :) = true; is_boundary(end, :) = true;
is_boundary(:, 1) = true; is_boundary(:, end) = true;
bnd_idx = find(is_boundary(:));

A_mod = A; A_mod(bnd_idx, :) = 0; A_mod(bnd_idx, bnd_idx) = speye(length(bnd_idx));
B_mod = B; B_mod(bnd_idx, :) = 0;
C_mod = C; C_mod(bnd_idx, :) = 0;

t = dt;
for n = 2:Nt-1
    t = t + dt;
    rhs = B_mod * uCurr + C_mod * uPrev;
    rhs(bnd_idx) = 0; 
    uNext = A_mod \ rhs;
    uPrev = uCurr;
    uCurr = uNext;
end

u_num = reshape(uCurr, Ny, Nx);
u_ref = u_exact(X, Y, t);
err = max(abs(u_num(:) - u_ref(:)));
fprintf("Max Error with damping at T=%.2f: %e\n", t, err);
