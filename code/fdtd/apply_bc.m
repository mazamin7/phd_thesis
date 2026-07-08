function [A_mod,B_mod,C_mod,rhsFun,leftTrace,rightTrace] = ...
    apply_bc(A,B,C,bcType,Nx,dx,dt,c,gamma,nu,q,r,gLeft,gRight)

rhsFun = @(t) zeros(Nx,1);

leftTrace  = @(u,v,t) NaN;
rightTrace = @(u,v,t) NaN;

%------------------------------------------------------------
% Copy interior matrices
%------------------------------------------------------------

A_mod = A;
B_mod = B;
C_mod = C;

%------------------------------------------------------------
% Boundary data
%------------------------------------------------------------

switch lower(bcType)

%% ============================================================
% DIRICHLET
% ============================================================

case 'dirichlet'

    A_mod(1,:) = 0;
    B_mod(1,:) = 0;
    C_mod(1,:) = 0;
    A_mod(1,1) = 1;

    A_mod(end,:) = 0;
    B_mod(end,:) = 0;
    C_mod(end,:) = 0;
    A_mod(end,end) = 1;

    rhsFun = @(t) dirichlet_rhs(t,Nx,dx,gLeft,gRight);

%% ============================================================
% NEUMANN (ghost)
% ============================================================

case 'neumannghost'

    A_mod(1,2) = 2*A_mod(1,2);
    B_mod(1,2) = 2*B_mod(1,2);
    C_mod(1,2) = 2*C_mod(1,2);

    A_mod(end,end-1) = 2*A_mod(end,end-1);
    B_mod(end,end-1) = 2*B_mod(end,end-1);
    C_mod(end,end-1) = 2*C_mod(end,end-1);

    rhsFun = @(t) neumannGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight);

%% ============================================================
% ROBIN (classic)
% ============================================================

case 'robinclassic'

    A_mod(1,:) = 0;
    B_mod(1,:) = 0;
    C_mod(1,:) = 0;
    
    A_mod(1,1) = q/dt + r/2;
    
    B_mod(1,1) = q/dt - 1/dx - r/2;
    B_mod(1,2) = 1/dx;

    A_mod(end,:) = 0;
    B_mod(end,:) = 0;
    C_mod(end,:) = 0;
    
    A_mod(end,end) = q/dt + r/2;
    
    B_mod(end,end)   = q/dt - 1/dx - r/2;
    B_mod(end,end-1) = 1/dx;

    rhsFun = @(t) robinClassic_rhs(t,Nx,dx,gLeft,gRight);

%% ============================================================
% ROBIN (BDF)
% ============================================================

case 'robinbdf'

    A_mod(1,:) = 0;
    B_mod(1,:) = 0;
    C_mod(1,:) = 0;

    A_mod(1,1) =  3/(2*dx)+3*q/(2*dt)+r;
    A_mod(1,2) = -4/(2*dx);
    A_mod(1,3) =  1/(2*dx);

    B_mod(1,1) =  4*q/(2*dt);
    C_mod(1,1) = -q/(2*dt);

    A_mod(end,:) = 0;
    B_mod(end,:) = 0;
    C_mod(end,:) = 0;

    A_mod(end,end)   =  3/(2*dx)+3*q/(2*dt)+r;
    A_mod(end,end-1) = -4/(2*dx);
    A_mod(end,end-2) =  1/(2*dx);

    B_mod(end,end) =  4*q/(2*dt);
    C_mod(end,end) = -q/(2*dt);

    rhsFun = @(t) robinBdf_rhs(t,Nx,dx,dt,gLeft,gRight);

%% ============================================================
% ROBIN (ghost)
% ============================================================

case 'robinghost'

    alpha = c^2*q*dt/dx;
    sigma2 = (c*dt/dx)^2;
    % Left Robin
    A_mod(1,2) = 2*A(1,2);
    B_mod(1,2) = 2*B(1,2);
    C_mod(1,2) = 2*C(1,2);

    A_mod(1,1) = A(1,1) + alpha;
    B_mod(1,1) = B(1,1) - 2*sigma2*r*dx;
    C_mod(1,1) = C(1,1) + alpha;

    % Right Robin
    A_mod(end,end-1) = 2*A(end,end-1);
    B_mod(end,end-1) = 2*B(end,end-1);
    C_mod(end,end-1) = 2*C(end,end-1);

    A_mod(end,end) = A(end,end) + alpha;
    B_mod(end,end) = B(end,end) - 2*sigma2*r*dx;
    C_mod(end,end) = C(end,end) + alpha;

    rhsFun = @(t) robinGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight);

    leftTrace = @(u,v,t) ...
        ((u(2) - ...
        (u(2) + 2*dx*(gLeft(t)-q*v(1)-r*u(1))))/(2*dx)) ...
        + q*v(1) + r*u(1);
    
    rightTrace = @(u,v,t) ...
        -(((u(end-1) + 2*dx*(gRight(t)-q*v(end)-r*u(end))) ...
        - u(end-1))/(2*dx)) ...
        + q*v(end) + r*u(end);

%% ============================================================
% ROBIN (ghost) - NEUMANN (ghost)
% ============================================================

case 'robinneumannghost'

    alpha  = c^2*q*dt/dx;
    sigma2 = (c*dt/dx)^2;

    % Left Robin
    A_mod(1,2) = 2*A(1,2);
    B_mod(1,2) = 2*B(1,2);
    C_mod(1,2) = 2*C(1,2);

    A_mod(1,1) = A(1,1) + alpha;
    B_mod(1,1) = B(1,1) - 2*sigma2*r*dx;
    C_mod(1,1) = C(1,1) + alpha;

    % Right Neumann
    A_mod(end,end-1) = 2*A_mod(end,end-1);
    B_mod(end,end-1) = 2*B_mod(end,end-1);
    C_mod(end,end-1) = 2*C_mod(end,end-1);

    rhsFun = @(t) robinNeumannGhost_rhs( ...
        c,t,Nx,dx,dt,nu,gLeft,gRight);
    
    leftTrace = @(u,v,t) ...
       ((u(2) - ...
        (u(2) + 2*dx*(gLeft(t)-q*v(1)-r*u(1))))/(2*dx)) ...
        + q*v(1) + r*u(1);
    
    rightTrace = @(u,v,t) NaN;

%% ============================================================
% NEUMANN (ghost) - ROBIN (ghost)
% ============================================================

case 'neumannrobinghost'

    alpha  = c^2*q*dt/dx;
    sigma2 = (c*dt/dx)^2;

    % Left Neumann
    A_mod(1,2) = 2*A_mod(1,2);
    B_mod(1,2) = 2*B_mod(1,2);
    C_mod(1,2) = 2*C_mod(1,2);

    % Right Robin
    A_mod(end,end-1) = 2*A(end,end-1);
    B_mod(end,end-1) = 2*B(end,end-1);
    C_mod(end,end-1) = 2*C_mod(end,end-1);

    A_mod(end,end) = A(end,end) + alpha;
    B_mod(end,end) = B(end,end) - 2*sigma2*r*dx;
    C_mod(end,end) = C(end,end) + alpha;

    rhsFun = @(t) neumannRobinGhost_rhs( ...
        c,t,Nx,dx,dt,nu,gLeft,gRight);

    leftTrace = @(u,v,t) NaN;
    
    rightTrace = @(u,v,t) ...
        -(((u(end-1) + 2*dx*(gRight(t)-q*v(end)-r*u(end))) ...
        - u(end-1))/(2*dx)) ...
        + q*v(end) + r*u(end);

%% ============================================================
% ROBIN (BDF) - NEUMANN (ghost)
% ============================================================

case 'robinbdfneumannghost'

    % Left Robin BDF
    A_mod(1,:) = 0; B_mod(1,:) = 0; C_mod(1,:) = 0;
    A_mod(1,1) =  3/(2*dx)+3*q/(2*dt)+r;
    A_mod(1,2) = -4/(2*dx);
    A_mod(1,3) =  1/(2*dx);
    B_mod(1,1) =  4*q/(2*dt);
    C_mod(1,1) = -q/(2*dt);

    % Right Neumann Ghost
    A_mod(end,end-1) = 2*A_mod(end,end-1);
    B_mod(end,end-1) = 2*B_mod(end,end-1);
    C_mod(end,end-1) = 2*C_mod(end,end-1);

    rhsFun = @(t) robinBdfNeumannGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight);
    
    leftTrace = @(u,v,t) NaN;
    rightTrace = @(u,v,t) NaN;

%% ============================================================
% NEUMANN (ghost) - ROBIN (BDF)
% ============================================================

case 'neumannghostrobinbdf'

    % Left Neumann Ghost
    A_mod(1,2) = 2*A_mod(1,2);
    B_mod(1,2) = 2*B_mod(1,2);
    C_mod(1,2) = 2*C_mod(1,2);

    % Right Robin BDF
    A_mod(end,:) = 0; B_mod(end,:) = 0; C_mod(end,:) = 0;
    A_mod(end,end)   =  3/(2*dx)+3*q/(2*dt)+r;
    A_mod(end,end-1) = -4/(2*dx);
    A_mod(end,end-2) =  1/(2*dx);
    B_mod(end,end) =  4*q/(2*dt);
    C_mod(end,end) = -q/(2*dt);

    rhsFun = @(t) neumannGhostRobinBdf_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight);

    leftTrace = @(u,v,t) NaN;
    rightTrace = @(u,v,t) NaN;

otherwise

    error('Unknown boundary condition.')

end

end

%% ============================================================
% RHS FUNCTIONS
% ============================================================

function rhs = dirichlet_rhs(t,Nx,dx,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1)   = gLeft(t);
rhs(end) = gRight(t);
end

function rhs = robinClassic_rhs(t,Nx,dx,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1)   = gLeft(t);
rhs(end) = gRight(t);
end

function rhs = robinBdf_rhs(t,Nx,dx,dt,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1)   = gLeft(t+dt);
rhs(end) = gRight(t+dt);
end

function rhs = neumannGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = -2*c^2*dt^2/dx * gLeft(t) - nu*dt/dx * (gLeft(t+dt) - gLeft(t-dt));
rhs(end) = 2*c^2*dt^2/dx * gRight(t) - nu*dt/dx * (gRight(t+dt) - gRight(t-dt));
end

function rhs = robinGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = 2*c^2*dt^2/dx * gLeft(t) - nu*dt/dx * (gLeft(t+dt) - gLeft(t-dt));
rhs(end) = 2*c^2*dt^2/dx * gRight(t) - nu*dt/dx * (gRight(t+dt) - gRight(t-dt));
end

function rhs = robinNeumannGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = 2*c^2*dt^2/dx * gLeft(t) - nu*dt/dx * (gLeft(t+dt) - gLeft(t-dt));
rhs(end) = 2*c^2*dt^2/dx * gRight(t) - nu*dt/dx * (gRight(t+dt) - gRight(t-dt));
end

function rhs = neumannRobinGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = -2*c^2*dt^2/dx * gLeft(t) - nu*dt/dx * (gLeft(t+dt) - gLeft(t-dt));
rhs(end) = 2*c^2*dt^2/dx * gRight(t) - nu*dt/dx * (gRight(t+dt) - gRight(t-dt));
end

function rhs = robinBdfNeumannGhost_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = gLeft(t+dt);
rhs(end) = 2*c^2*dt^2/dx * gRight(t) - nu*dt/dx * (gRight(t+dt) - gRight(t-dt));
end

function rhs = neumannGhostRobinBdf_rhs(c,t,Nx,dx,dt,nu,gLeft,gRight)
rhs = zeros(Nx,1);
rhs(1) = -2*c^2*dt^2/dx * gLeft(t) - nu*dt/dx * (gLeft(t+dt) - gLeft(t-dt));
rhs(end) = gRight(t+dt);
end
