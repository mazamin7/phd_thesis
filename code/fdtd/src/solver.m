function [x_grid,t_grid,u_array,v_array,robinLeft,robinRight] = solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    bcType,q,r,...
    gLeft,gRight)

%-------------------------------------------------------------
% Default boundary waveforms
%-------------------------------------------------------------

if nargin < 14
    gLeft = @(t) 0;
    gRight = @(t) 0;
end

%-------------------------------------------------------------
% Grid
%-------------------------------------------------------------

x_grid = (0:dx:L).';
t_grid = (0:dt:T);

Nx = length(x_grid);
Nt = length(t_grid);

u_array = zeros(Nx,Nt);
v_array = zeros(Nx,Nt);

%-------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------

u0 = u0_fun(x_grid);
v0 = v0_fun(x_grid);

u_array(:,1) = u0;
v_array(:,1) = v0;

%-------------------------------------------------------------
% Robin traces
%-------------------------------------------------------------

robinLeft  = zeros(1,Nt);
robinRight = zeros(1,Nt);

%-------------------------------------------------------------
% Taylor initialization
%-------------------------------------------------------------

uxx0 = zeros(size(x_grid));

for i = 2:Nx-1
    uxx0(i) = (u0(i+1)-2*u0(i)+u0(i-1))/dx^2;
end

u0_left = u0_fun(x_grid(1) - dx);
u0_right = u0_fun(x_grid(end) + dx);
uxx0(1)   = (u0(2) - 2*u0(1) + u0_left)/dx^2;
uxx0(end) = (u0_right - 2*u0(end) + u0(end-1))/dx^2;

vxx0 = zeros(size(x_grid));

for i = 2:Nx-1
    vxx0(i) = (v0(i+1)-2*v0(i)+v0(i-1))/dx^2;
end

v0_left = v0_fun(x_grid(1) - dx);
v0_right = v0_fun(x_grid(end) + dx);
vxx0(1)   = (v0(2) - 2*v0(1) + v0_left)/dx^2;
vxx0(end) = (v0_right - 2*v0(end) + v0(end-1))/dx^2;

utt0 = ...
      c^2*uxx0 ...
    + nu*vxx0 ...
    - gamma*v0 ...
    + f_fun(x_grid,0);

uPrev = u0;

uCurr = ...
      u0 ...
    + dt*v0 ...
    + 0.5*dt^2*utt0;

u_array(:,2) = uCurr;

% second-order approximation of the velocity at t = dt
v_array(:,2) = v0 + dt*utt0;

%-------------------------------------------------------------
% Interior matrices
%-------------------------------------------------------------

[A,B,C] = interior_assembly( ...
    Nx,dx,dt,...
    c,gamma,nu);

%-------------------------------------------------------------
% Boundary conditions
%-------------------------------------------------------------

[A_mod,B_mod,C_mod,rhsFunBc,leftTrace,rightTrace] = apply_bc( ...
    A,B,C,...
    bcType,Nx,dx,dt,...
    c,gamma,nu,...
    q,r,...
    gLeft,...
    gRight);

%-------------------------------------------------------------
% Time stepping
%-------------------------------------------------------------

for n = 2:Nt-1

    t = t_grid(n);

    rhs = ...
          B_mod*uCurr ...
        + C_mod*uPrev ...
        + dt^2*f_fun(x_grid,t);

    rhs = rhs + rhsFunBc(t);

    uNext = A_mod\rhs;

    u_array(:,n+1) = uNext;

    % centered approximation of velocity
    v_array(:,n) = (uNext-uPrev)/(2*dt);

    uPrev = uCurr;
    uCurr = uNext;

end

%-------------------------------------------------------------
% Last velocity
%-------------------------------------------------------------

v_array(:,Nt) = ...
    (3*u_array(:,Nt)-4*u_array(:,Nt-1)+u_array(:,Nt-2))/(2*dt);

%-------------------------------------------------------------
% Robin traces
%-------------------------------------------------------------

for n = 1:Nt

    robinLeft(n) = leftTrace( ...
        u_array(:,n), ...
        v_array(:,n), ...
        t_grid(n));

    robinRight(n) = rightTrace( ...
        u_array(:,n), ...
        v_array(:,n), ...
        t_grid(n));

end

end