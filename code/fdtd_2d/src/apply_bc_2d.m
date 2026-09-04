function [A_mod,B_mod,C_mod,rhsFunBc] = apply_bc_2d( ...
    A, B, C, Nx, Ny, dx, dy, dt, q, r, ...
    bcLeft, bcRight, bcBottom, bcTop, ...
    gLeft, gRight, gBottom, gTop)

    A_mod = A; B_mod = B; C_mod = C;
    
    % Function to build RHS
    rhsFunBc = @(t) build_rhs(t, Nx, Ny, dx, dt, bcLeft, bcRight, bcBottom, bcTop, gLeft, gRight, gBottom, gTop);

    % Indices of boundaries
    % y is fast (1:Ny), x is slow (1:Nx)
    idxLeft   = 1 : Ny;
    idxRight  = (Nx-1)*Ny + 1 : Nx*Ny;
    idxBottom = 1 : Ny : Nx*Ny; % y=0
    idxTop    = Ny : Ny : Nx*Ny; % y=Ly
    
    % --- LEFT BOUNDARY ---
    if strcmp(bcLeft, 'dirichlet')
        A_mod(idxLeft, :) = 0; B_mod(idxLeft, :) = 0; C_mod(idxLeft, :) = 0;
        A_mod(idxLeft, idxLeft) = speye(Ny);
    elseif strcmp(bcLeft, 'robinbdf')
        A_mod(idxLeft, :) = 0; B_mod(idxLeft, :) = 0; C_mod(idxLeft, :) = 0;
        for k = idxLeft
            A_mod(k, k) = 3/(2*dx) + 3*q/(2*dt) + r;
            A_mod(k, k+Ny) = -4/(2*dx);
            A_mod(k, k+2*Ny) = 1/(2*dx);
            B_mod(k, k) = 4*q/(2*dt);
            C_mod(k, k) = -q/(2*dt);
        end
    end
    
    % --- RIGHT BOUNDARY ---
    if strcmp(bcRight, 'dirichlet')
        A_mod(idxRight, :) = 0; B_mod(idxRight, :) = 0; C_mod(idxRight, :) = 0;
        A_mod(idxRight, idxRight) = speye(Ny);
    elseif strcmp(bcRight, 'robinbdf')
        A_mod(idxRight, :) = 0; B_mod(idxRight, :) = 0; C_mod(idxRight, :) = 0;
        for k = idxRight
            A_mod(k, k) = 3/(2*dx) + 3*q/(2*dt) + r;
            A_mod(k, k-Ny) = -4/(2*dx);
            A_mod(k, k-2*Ny) = 1/(2*dx);
            B_mod(k, k) = 4*q/(2*dt);
            C_mod(k, k) = -q/(2*dt);
        end
    end
    
    % --- BOTTOM BOUNDARY --- (y=0)
    if strcmp(bcBottom, 'dirichlet')
        A_mod(idxBottom, :) = 0; B_mod(idxBottom, :) = 0; C_mod(idxBottom, :) = 0;
        A_mod(idxBottom, idxBottom) = speye(Nx);
    end
    
    % --- TOP BOUNDARY --- (y=Ly)
    if strcmp(bcTop, 'dirichlet')
        A_mod(idxTop, :) = 0; B_mod(idxTop, :) = 0; C_mod(idxTop, :) = 0;
        A_mod(idxTop, idxTop) = speye(Nx);
    end
    
end

function rhs = build_rhs(t, Nx, Ny, dx, dt, bcLeft, bcRight, bcBottom, bcTop, gLeft, gRight, gBottom, gTop)
    rhs = zeros(Nx * Ny, 1);
    
    idxLeft   = 1 : Ny;
    idxRight  = (Nx-1)*Ny + 1 : Nx*Ny;
    idxBottom = 1 : Ny : Nx*Ny;
    idxTop    = Ny : Ny : Nx*Ny;
    
    if strcmp(bcLeft, 'dirichlet')
        rhs(idxLeft) = gLeft(t);
    elseif strcmp(bcLeft, 'robinbdf')
        rhs(idxLeft) = gLeft(t+dt);
    end
    
    if strcmp(bcRight, 'dirichlet')
        rhs(idxRight) = gRight(t);
    elseif strcmp(bcRight, 'robinbdf')
        rhs(idxRight) = gRight(t+dt);
    end
    
    if strcmp(bcBottom, 'dirichlet')
        rhs(idxBottom) = gBottom(t);
    end
    
    if strcmp(bcTop, 'dirichlet')
        rhs(idxTop) = gTop(t);
    end
end
