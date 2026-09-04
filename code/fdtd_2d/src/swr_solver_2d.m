function [X, Y, t_grid, u_iter, v_iter] = swr_solver_2d( ...
    u0_fun, v0_fun, f_fun, ...
    dx, dy, dt, Lx, Ly, T_end, ...
    c, gamma, nu, q, r, ...
    delta, iterations)

    x_grid = (0:dx:Lx)';
    y_grid = (0:dy:Ly)';
    t_grid = 0:dt:T_end;

    Nx = length(x_grid);
    Ny = length(y_grid);
    Nt = length(t_grid);
    [X, Y] = meshgrid(x_grid, y_grid);

    u_iter = zeros(Ny, Nx, Nt, iterations);
    v_iter = zeros(Ny, Nx, Nt, iterations);

    % ALIGN MID WITH GRID
    mid = round((Lx/2)/dx) * dx;
    delta_half_L = ceil( (delta/2)/dx ) * dx;
    delta_half_R = floor( (delta/2)/dx ) * dx;

    xL0 = 0; xL1 = mid + delta_half_L;
    xR0 = mid - delta_half_R; xR1 = Lx;

    xL = (xL0:dx:xL1)'; xR = (xR0:dx:xR1)';
    NL = length(xL); NR = length(xR);

    gLL = @(t) zeros(Ny, 1);
    gLR = @(t) zeros(Ny, 1);
    gRL = @(t) zeros(Ny, 1);
    gRR = @(t) zeros(Ny, 1);

    for k = 1:iterations
        % Left
        [~, ~, ~, uL, vL] = solver_2d( ...
            @(X_sub, Y_sub) u0_fun(X_sub, Y_sub), ...
            @(X_sub, Y_sub) v0_fun(X_sub, Y_sub), ...
            @(X_sub, Y_sub, t) f_fun(X_sub, Y_sub, t), ...
            dx, dy, dt, xL1-xL0, Ly, T_end, ...
            c, gamma, nu, q, r, ...
            'dirichlet', 'robinbdf', 'dirichlet', 'dirichlet', ...
            gLL, gRL, @(t) zeros(NL,1), @(t) zeros(NL,1));

        % Right
        [~, ~, ~, uR, vR] = solver_2d( ...
            @(X_sub, Y_sub) u0_fun(X_sub + xR0, Y_sub), ...
            @(X_sub, Y_sub) v0_fun(X_sub + xR0, Y_sub), ...
            @(X_sub, Y_sub, t) f_fun(X_sub + xR0, Y_sub, t), ...
            dx, dy, dt, xR1-xR0, Ly, T_end, ...
            c, gamma, nu, q, r, ...
            'robinbdf', 'dirichlet', 'dirichlet', 'dirichlet', ...
            gLR, gRR, @(t) zeros(NR,1), @(t) zeros(NR,1));

        u = zeros(Ny, Nx, Nt);
        v = zeros(Ny, Nx, Nt);
        
        for i = 1:Nx
            if x_grid(i) <= mid
                j = round((x_grid(i)-xL0)/dx)+1;
                j = max(1, min(j, NL));
                u(:, i, :) = uL(:, j, :);
                v(:, i, :) = vL(:, j, :);
            else
                j = round((x_grid(i)-xR0)/dx)+1;
                j = max(1, min(j, NR));
                u(:, i, :) = uR(:, j, :);
                v(:, i, :) = vR(:, j, :);
            end
        end

        u_iter(:,:,:,k) = u;
        v_iter(:,:,:,k) = v;

        idxL = round((xR0 - xL0)/dx) + 1;
        idxL = max(1, min(idxL, NL));
        robinLeftVals = zeros(Ny, Nt);
        for n = 3:Nt
            if round(delta/dx) == 1
                uxL_neg = (uL(:, idxL, n) - uL(:, idxL+1, n)) / dx;
            else
                uxL_neg = (3*uL(:, idxL, n) - 4*uL(:, idxL+1, n) + uL(:, idxL+2, n)) / (2*dx);
            end
            utL = (3*uL(:, idxL, n) - 4*uL(:, idxL, n-1) + uL(:, idxL, n-2)) / (2*dt);
            robinLeftVals(:, n) = uxL_neg + q*utL + r*uL(:, idxL, n);
        end
        robinLeftVals(:, 1:2) = repmat(robinLeftVals(:, 3), 1, 2);

        idxR = round((xL1 - xR0)/dx) + 1;
        idxR = max(1, min(idxR, NR));
        robinRightVals = zeros(Ny, Nt);
        for n = 3:Nt
            if round(delta/dx) == 1
                uxR = (uR(:, idxR, n) - uR(:, idxR-1, n)) / dx;
            else
                uxR = (3*uR(:, idxR, n) - 4*uR(:, idxR-1, n) + uR(:, idxR-2, n)) / (2*dx);
            end
            utR = (3*uR(:, idxR, n) - 4*uR(:, idxR, n-1) + uR(:, idxR, n-2)) / (2*dt);
            robinRightVals(:, n) = uxR + q*utR + r*uR(:, idxR, n);
        end
        robinRightVals(:, 1:2) = repmat(robinRightVals(:, 3), 1, 2);

        gLR = @(t) interp1(t_grid, robinLeftVals', t, 'linear', 'extrap')';
        gRL = @(t) interp1(t_grid, robinRightVals', t, 'linear', 'extrap')';
    end
end
