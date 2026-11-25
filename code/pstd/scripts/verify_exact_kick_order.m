% Verify temporal accuracy of the exact-kick ARD solver against the exact
% solution of the same semidiscrete ARD operator.
clear; clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));

L = 1;
c = 1;
gamma = 0.03;
nu = 0.05;
N = 32;
dx = L/N;
T = 0.1;
space_orders = [2 4];
step_counts = [10 20 40 80];
dt_values = T./step_counts;

bcType = 'neumannGhost'; %#ok<NASGU>
q = 1/c; %#ok<NASGU>
r = 0; %#ok<NASGU>
u0_fun = @(x) exp(-90*(x-0.37).^2) + 0.2*cos(3*pi*x/L);
v0_fun = @(x) 0.15*sin(2*pi*x/L);
f_fun = @(x,t) zeros(size(x)); %#ok<NASGU>

NxSub = N/2;
xL = ((0:NxSub-1)'+0.5)*dx;
xR = L/2 + ((0:NxSub-1)'+0.5)*dx;
x = [xL; xR];
u0 = u0_fun(x);
v0 = v0_fun(x);

% Build the exact matrix of the semidiscrete equation followed by both ARD
% solvers. The block DCT operator is obtained by applying it to basis vectors.
modes = (0:NxSub-1)';
lambda = -(modes*pi/(L/2)).^2;
Q = dct(eye(NxSub));
Llocal = idct(diag(lambda)*Q);
Ldrift = blkdiag(Llocal,Llocal);

for space_order = space_orders
    C_res = get_residue_matrix(N,space_order);
    Ltotal = Ldrift + C_res/dx^2;
    A = [zeros(N), eye(N); ...
         c^2*Ltotal, -gamma*eye(N) + nu*Ltotal];
    y_ref = expm(T*A)*[u0;v0];
    u_ref = y_ref(1:N);

    err_euler = zeros(size(dt_values));
    err_exact = zeros(size(dt_values));
    for k = 1:numel(dt_values)
        dt = dt_values(k);
        [~,~,~,~,u_euler] = ard_solver( ...
            u0_fun,v0_fun,f_fun,dx,dt,L,T,c,gamma,nu, ...
            'neumannGhost',1/c,0,space_order);
        [~,~,~,~,u_exact] = ard_solver_exact_kick( ...
            u0_fun,v0_fun,f_fun,dx,dt,L,T,c,gamma,nu, ...
            'neumannGhost',1/c,0,space_order);

        err_euler(k) = norm(u_euler(:,end)-u_ref,inf);
        err_exact(k) = norm(u_exact(:,end)-u_ref,inf);
    end

    rate_euler = log2(err_euler(1:end-1)./err_euler(2:end));
    rate_exact = log2(err_exact(1:end-1)./err_exact(2:end));
    fprintf('\nInterface order %d, N=%d, T=%.3g\n',space_order,N,T);
    fprintf('       dt          Euler error    Euler rate     Exact-kick error  Exact-kick rate\n');
    for k = 1:numel(dt_values)
        if k < numel(dt_values)
            fprintf('%11.4g   %13.5e   %10.4f   %16.5e   %14.4f\n', ...
                dt_values(k),err_euler(k),rate_euler(k),err_exact(k),rate_exact(k));
        else
            fprintf('%11.4g   %13.5e        --     %16.5e        --\n', ...
                dt_values(k),err_euler(k),err_exact(k));
        end
    end
end
