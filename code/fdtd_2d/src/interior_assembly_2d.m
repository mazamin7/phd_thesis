function [A,B,C] = interior_assembly_2d(Nx, Ny, dx, dy, dt, c, gamma, nu)
    % Identity
    Ix = speye(Nx);
    Iy = speye(Ny);
    I = speye(Nx * Ny);

    % Discrete Laplacian
    e_x = ones(Nx,1);
    e_y = ones(Ny,1);

    Dxx = spdiags([e_x -2*e_x e_x],[-1 0 1],Nx,Nx) / dx^2;
    Dyy = spdiags([e_y -2*e_y e_y],[-1 0 1],Ny,Ny) / dy^2;

    % 2D Laplacian using Kronecker tensor product
    % Assuming flattening stacks columns (y is fast index, x is slow index)
    D = kron(Dxx, Iy) + kron(Ix, Dyy);

    % Matrices
    A = (1 + gamma*dt/2)*I - (nu*dt/2)*D;
    B = 2*I + c^2*dt^2*D;
    C = -(1 - gamma*dt/2)*I - (nu*dt/2)*D;
end
