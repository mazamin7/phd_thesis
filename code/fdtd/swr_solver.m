function [x_grid,t_grid,u_iter,v_iter] = swr_solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    q,r,...
    delta,...
    iterations)

%-------------------------------------------------------------
% Global grid
%-------------------------------------------------------------

x_grid = (0:dx:L).';
t_grid = (0:dt:T);

Nx = length(x_grid);
Nt = length(t_grid);

u_iter = zeros(Nx,Nt,iterations);
v_iter = zeros(Nx,Nt,iterations);

%-------------------------------------------------------------
% Overlapping decomposition
%-------------------------------------------------------------

mid = L/2;

delta_half_L = ceil( (delta/2)/dx ) * dx;
delta_half_R = floor( (delta/2)/dx ) * dx;

xL0 = 0;
xL1 = mid + delta_half_L;

xR0 = mid - delta_half_R;
xR1 = L;

xL = (xL0:dx:xL1).';
xR = (xR0:dx:xR1).';

NL = length(xL);
NR = length(xR);

% indices in the global mesh

iL = 1:NL;

firstRight = round(xR0/dx)+1;
iR = firstRight:Nx;

% midpoint used for assembling

midIndex = round(mid/dx)+1;

%-------------------------------------------------------------
% Initial transmission data
%-------------------------------------------------------------

gLL = @(t) 0;
gLR = @(t) 0;
gRL = @(t) 0;
gRR = @(t) 0;

%-------------------------------------------------------------
% Schwarz iterations
%-------------------------------------------------------------

for k = 1:iterations

    fprintf('Iteration %d\n',k);

    %=========================================================
    % Left subdomain
    %=========================================================

    [~,~,uL,vL,~,robinLeft] = solver( ...
        @(x) u0_fun(x),...
        @(x) v0_fun(x),...
        @(x,t) f_fun(x,t),...
        dx,dt,xL1-xL0,T,...
        c,gamma,nu,...
        'neumannGhostRobinBdf',...
        q,r,...
        gLL,...
        gRL);

    %=========================================================
    % Right subdomain
    %=========================================================

    [~,~,uR,vR,robinRight,~] = solver( ...
        @(x) u0_fun(x+xR0),...
        @(x) v0_fun(x+xR0),...
        @(x,t) f_fun(x+xR0,t),...
        dx,dt,xR1-xR0,T,...
        c,gamma,nu,...
        'robinBdfNeumannGhost',...
        q,r,...
        gLR,...
        gRR);

    %=========================================================
    % Assemble global solution
    %=========================================================

    u = zeros(Nx,Nt);
    v = zeros(Nx,Nt);
    
    for i = 1:Nx
    
        if x_grid(i) <= mid
    
            j = round((x_grid(i)-xL0)/dx)+1;
    
            u(i,:) = uL(j,:);
            v(i,:) = vL(j,:);
    
        else
    
            j = round((x_grid(i)-xR0)/dx)+1;
    
            u(i,:) = uR(j,:);
            v(i,:) = vR(j,:);
    
        end
    
    end


    u_iter(:,:,k) = u;
    v_iter(:,:,k) = v;

    %=========================================================
    % Exchange Robin traces
    %=========================================================

    if delta > 0
        % Left subdomain sends data to Right subdomain's left boundary (xR0)
        idxL = round((xR0 - xL0)/dx) + 1;
        robinLeftVals = zeros(1, Nt);
        for n = 3:Nt
            if round(delta/dx) == 1
                uxL_neg = (uL(idxL,n) - uL(idxL+1,n))/dx;
            else
                uxL_neg = (3*uL(idxL,n) - 4*uL(idxL+1,n) + uL(idxL+2,n))/(2*dx);
            end
            utL = (3*uL(idxL,n) - 4*uL(idxL,n-1) + uL(idxL,n-2))/(2*dt);
            robinLeftVals(n) = uxL_neg + q*utL + r*uL(idxL,n);
        end
        robinLeftVals(1:2) = robinLeftVals(3);

        % Right subdomain sends data to Left subdomain's right boundary (xL1)
        idxR = round((xL1 - xR0)/dx) + 1;
        robinRightVals = zeros(1, Nt);
        for n = 3:Nt
            if round(delta/dx) == 1
                uxR = (uR(idxR,n) - uR(idxR-1,n))/dx;
            else
                uxR = (3*uR(idxR,n) - 4*uR(idxR-1,n) + uR(idxR-2,n))/(2*dx);
            end
            utR = (3*uR(idxR,n) - 4*uR(idxR,n-1) + uR(idxR,n-2))/(2*dt);
            robinRightVals(n) = uxR + q*utR + r*uR(idxR,n);
        end
        robinRightVals(1:2) = robinRightVals(3);
    else
        robinLeftVals  = robinLeft;
        robinRightVals = robinRight;
    end

    gLR = @(t) interp1( ...
        t_grid,...
        robinLeftVals,...
        t,...
        'linear',...
        'extrap');

    gRL = @(t) interp1( ...
        t_grid,...
        robinRightVals,...
        t,...
        'linear',...
        'extrap');

end

disp('finished SWR iterations')

end