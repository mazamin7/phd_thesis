function [X, Y, t_grid, u_array, v_array, robinLeft, robinRight] = solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, q, r, ...
    bcLeft, bcRight, bcBottom, bcTop, ...
    gLeft, gRight, gBottom, gTop)

    x = (0:dx:Lx)';
    y = (0:dy:Ly)';
    t_grid = 0:dt:T_end;
    Nt = length(t_grid);

    Nx = length(x); Ny = length(y);
    [X, Y] = meshgrid(x, y);

    u_array = zeros(Ny, Nx, Nt);
    v_array = zeros(Ny, Nx, Nt);
    robinLeft = zeros(Ny, Nt);
    robinRight = zeros(Ny, Nt);

    u0 = u0_fun(X, Y); v0 = v0_fun(X, Y);
    u_array(:, :, 1) = u0; v_array(:, :, 1) = v0;
    
    u0_left = u0_fun(X(:,1) - dx, Y(:,1));
    u0_right = u0_fun(X(:,end) + dx, Y(:,end));
    u0_bottom = u0_fun(X(1,:), Y(1,:) - dy);
    u0_top = u0_fun(X(end,:), Y(end,:) + dy);
    
    v0_left = v0_fun(X(:,1) - dx, Y(:,1));
    v0_right = v0_fun(X(:,end) + dx, Y(:,end));
    v0_bottom = v0_fun(X(1,:), Y(1,:) - dy);
    v0_top = v0_fun(X(end,:), Y(end,:) + dy);
    
    uxx0 = zeros(Ny, Nx); uyy0 = zeros(Ny, Nx);
    vxx0 = zeros(Ny, Nx); vyy0 = zeros(Ny, Nx);
    
    for i = 2:Nx-1
        uxx0(:, i) = (u0(:, i+1) - 2*u0(:, i) + u0(:, i-1)) / dx^2;
        vxx0(:, i) = (v0(:, i+1) - 2*v0(:, i) + v0(:, i-1)) / dx^2;
    end
    uxx0(:, 1) = (u0(:, 2) - 2*u0(:, 1) + u0_left) / dx^2;
    uxx0(:, end) = (u0_right - 2*u0(:, end) + u0(:, end-1)) / dx^2;
    
    vxx0(:, 1) = (v0(:, 2) - 2*v0(:, 1) + v0_left) / dx^2;
    vxx0(:, end) = (v0_right - 2*v0(:, end) + v0(:, end-1)) / dx^2;
    
    for j = 2:Ny-1
        uyy0(j, :) = (u0(j+1, :) - 2*u0(j, :) + u0(j-1, :)) / dy^2;
        vyy0(j, :) = (v0(j+1, :) - 2*v0(j, :) + v0(j-1, :)) / dy^2;
    end
    uyy0(1, :) = (u0(2, :) - 2*u0(1, :) + u0_bottom) / dy^2;
    uyy0(end, :) = (u0_top - 2*u0(end, :) + u0(end-1, :)) / dy^2;
    
    vyy0(1, :) = (v0(2, :) - 2*v0(1, :) + v0_bottom) / dy^2;
    vyy0(end, :) = (v0_top - 2*v0(end, :) + v0(end-1, :)) / dy^2;
    
    Du0 = uxx0 + uyy0;
    Dv0 = vxx0 + vyy0;
    
    f0 = f_fun(X, Y, 0);
    utt0 = c^2 * Du0 + nu * Dv0 - gamma * v0 + f0;

    uCurr = u0(:);
    uNext = uCurr + dt * v0(:) + 0.5 * dt^2 * utt0(:);
    vNext = v0(:) + dt * utt0(:);

    [A, B, C] = interior_assembly_2d(Nx, Ny, dx, dy, dt, c, gamma, nu);

    [A_mod, B_mod, C_mod, rhsFunBc] = apply_bc_2d( ...
        A, B, C, Nx, Ny, dx, dy, dt, q, r, ...
        bcLeft, bcRight, bcBottom, bcTop, ...
        gLeft, gRight, gBottom, gTop);

    idxLeft   = 1 : Ny;
    idxRight  = (Nx-1)*Ny + 1 : Nx*Ny;
    idxBottom = 1 : Ny : Nx*Ny;
    idxTop    = Ny : Ny : Nx*Ny;
    
    if strcmp(bcLeft, 'dirichlet'), uNext(idxLeft) = gLeft(dt); vNext(idxLeft) = (gLeft(dt)-gLeft(0))/dt; end
    if strcmp(bcRight, 'dirichlet'), uNext(idxRight) = gRight(dt); vNext(idxRight) = (gRight(dt)-gRight(0))/dt; end
    if strcmp(bcBottom, 'dirichlet'), uNext(idxBottom) = gBottom(dt); vNext(idxBottom) = (gBottom(dt)-gBottom(0))/dt; end
    if strcmp(bcTop, 'dirichlet'), uNext(idxTop) = gTop(dt); vNext(idxTop) = (gTop(dt)-gTop(0))/dt; end

    u_array(:, :, 2) = reshape(uNext, Ny, Nx);
    v_array(:, :, 2) = reshape(vNext, Ny, Nx);

    uPrev = uCurr;
    uCurr = uNext;

    for n = 2:Nt-1
        t = t_grid(n);
        f_curr = reshape(f_fun(X, Y, t), [], 1);
        rhs = B_mod * uCurr + C_mod * uPrev + dt^2 * f_curr;
        rhs = rhs + rhsFunBc(t);

        uNext = A_mod \ rhs;
        u_array(:, :, n+1) = reshape(uNext, Ny, Nx);
        v_array(:, :, n) = reshape((uNext - uPrev) / (2 * dt), Ny, Nx);

        uPrev = uCurr;
        uCurr = uNext;
    end

    u_last = u_array(:,:,Nt); u_prev = u_array(:,:,Nt-1); u_prev2 = u_array(:,:,Nt-2);
    v_array(:, :, Nt) = (3*u_last - 4*u_prev + u_prev2) / (2 * dt);
end
