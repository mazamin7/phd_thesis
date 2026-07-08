clear all; close all; clc;
q = 1;
r = 0;
iterations = 10;
delta = 0.1;
L = 1;
T = 2;
dx = 0.01;
dt = 0.01;
c = 1;
gamma = 0;
nu = 0;
k1 = pi/L;
uExact = @(x,t) cos(k1*x).*cos(k1*c*t);
u0_fun = @(x) uExact(x,0);
v0_fun = @(x) zeros(size(x));
f_fun = @(x,t) zeros(size(x));
[x_grid,t_grid,u_iter,v_iter] = swr_solver(u0_fun,v0_fun,f_fun,dx,dt,L,T,c,gamma,nu,q,r,delta,iterations);
err = zeros(iterations,1);
for k=1:iterations
    err(k) = norm(u_iter(:,end,k) - uExact(x_grid, T));
end
disp('Errors:');
disp(err);
