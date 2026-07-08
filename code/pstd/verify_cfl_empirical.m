% verify_cfl_empirical.m
clear all; close all; clc;
addpath('WASAbi1D/utils');

orders = [2, 4, 6, 8];
computed_limits = [1.1708, 1.0599, 1.0180, 0.9918];
epsilon = 0.05;

L = 1;
Nx = 100;
dx = L/Nx;
c = 1;
T_sim = 250.0; % Very long time to let mild instabilities grow
gamma = 0;
nu = 0;
bcType = 'neumann';
q = 0; r = 0;

u0_fun = @(x) exp(-400*(x-L/2).^2);
v0_fun = @(x) zeros(size(x));
f_fun = @(x,t) zeros(size(x));

results = struct();

for k = 1:length(orders)
    space_order = orders(k);
    limit = computed_limits(k);
    
    cfl_stable = limit - epsilon;
    cfl_unstable = limit + epsilon;
    
    disp(['Testing Order ', num2str(space_order), '...']);
    
    % --- Test Stable ---
    dt_stable = cfl_stable * dx / c;
    try
        [~,~,~,~,u_stable,~,~,~,~] = ard_solver(u0_fun,v0_fun,f_fun,dx,dt_stable,L,T_sim,c,gamma,nu,bcType,q,r,space_order);
        max_val_stable = max(abs(u_stable(:)));
    catch ME
        disp(['Stable run error: ', ME.message]);
        max_val_stable = Inf;
    end
    
    % --- Test Unstable ---
    dt_unstable = cfl_unstable * dx / c;
    try
        [~,~,~,~,u_unstable,~,~,~,~] = ard_solver(u0_fun,v0_fun,f_fun,dx,dt_unstable,L,T_sim,c,gamma,nu,bcType,q,r,space_order);
        max_val_unstable = max(abs(u_unstable(:)));
    catch ME
        disp(['Unstable run error: ', ME.message]);
        max_val_unstable = Inf;
    end
    
    disp(['  CFL = ', num2str(cfl_stable), ' (Stable Expected): Max Amplitude = ', num2str(max_val_stable, '%e')]);
    disp(['  CFL = ', num2str(cfl_unstable), ' (Unstable Expected): Max Amplitude = ', num2str(max_val_unstable, '%e')]);
    
    results(k).order = space_order;
    results(k).stable_amp = max_val_stable;
    results(k).unstable_amp = max_val_unstable;
end

disp('Done!');
