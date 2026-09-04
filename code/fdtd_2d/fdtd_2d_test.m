Lx = pi;
Ly = pi;
Nx = 50;
Ny = 50;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);

c = 1.0;
gamma = 0.0;
nu = 0.0;

dt = 0.5 * min(dx, dy) / c;
T_end = 2.0;
Nt = floor(T_end / dt) + 1;
dt = T_end / (Nt - 1);

x = linspace(0, Lx, Nx)';
y = linspace(0, Ly, Ny)';
[X, Y] = meshgrid(x, y);

% Analytical solution:
% u(x,y,t) = cos(sqrt(2)*c*t) * sin(x) * sin(y)
u_exact = @(X, Y, t) cos(sqrt(2)*c*t) .* sin(X) .* sin(Y);
v_exact = @(X, Y, t) -sqrt(2)*c*sin(sqrt(2)*c*t) .* sin(X) .* sin(Y);
utt_exact = @(X, Y, t) -2*c^2*cos(sqrt(2)*c*t) .* sin(X) .* sin(Y);

% Flatten grid
u_flat = @(t) reshape(u_exact(X, Y, t), [], 1);

u0 = u_exact(X, Y, 0);
v0 = v_exact(X, Y, 0);
utt0 = utt_exact(X, Y, 0);

uCurr = u0(:);
uPrev = uCurr; % For n=1, we construct uNext = uCurr + dt*v0 + 0.5*dt^2*utt0
uNext = uCurr + dt*v0(:) + 0.5*dt^2*utt0(:);

uPrev = uCurr;
uCurr = uNext;

% Matrices
Ix = speye(Nx);
Iy = speye(Ny);
e_x = ones(Nx, 1);
e_y = ones(Ny, 1);
Dxx = spdiags([e_x -2*e_x e_x], [-1 0 1], Nx, Nx) / dx^2;
Dyy = spdiags([e_y -2*e_y e_y], [-1 0 1], Ny, Ny) / dy^2;

% Boundary conditions (Dirichlet)
% Zero out the boundary rows in Dxx and Dyy, but actually we can just apply BCs after.
% To make it easier, we assemble full D and apply BC to A, B, C.
D = kron(Dxx, Iy) + kron(Ix, Dyy);
I = speye(Nx * Ny);

A = (1 + gamma*dt/2)*I - (nu*dt/2)*D;
B = 2*I + c^2*dt^2*D;
C = -(1 - gamma*dt/2)*I - (nu*dt/2)*D;

% Identify boundary nodes
is_boundary = false(Ny, Nx);
is_boundary(1, :) = true;
is_boundary(end, :) = true;
is_boundary(:, 1) = true;
is_boundary(:, end) = true;
bnd_idx = find(is_boundary(:));
int_idx = find(~is_boundary(:));

% Apply Dirichlet BC to matrices
A_mod = A;
A_mod(bnd_idx, :) = 0;
A_mod(bnd_idx, bnd_idx) = speye(length(bnd_idx));

B_mod = B;
B_mod(bnd_idx, :) = 0;

C_mod = C;
C_mod(bnd_idx, :) = 0;

% Time loop
t = dt;
for n = 2:Nt-1
    t = t + dt;
    
    rhs = B_mod * uCurr + C_mod * uPrev;
    % Dirichlet boundary conditions are 0, so RHS at boundaries is 0
    rhs(bnd_idx) = 0; 
    
    uNext = A_mod \ rhs;
    
    uPrev = uCurr;
    uCurr = uNext;
end

u_num = reshape(uCurr, Ny, Nx);
u_ref = u_exact(X, Y, t);

err = max(abs(u_num(:) - u_ref(:)));
fprintf("Max Error at T=%.2f: %e\n", t, err);

