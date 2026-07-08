addpath(genpath('../src'));
% clear all; close all; clc;

%% ================================================================
% Sweep over Kelvin--Voigt viscosity nu: compare classical vs
% optimized Robin parameters, for a fixed experiment.
% ================================================================

%% Fixed parameters
c   = 1;
dx  = 0.02;
dt  = 0.02;

experiment = 'trapezoidStandingWaveNeumann';
bcType     = 'neumannGhost';
L = 1;
T = 2;
gamma = 0.5;

iterations = 30;
delta = 2*dx;
assert(delta >= 2*dx)
a = (L-delta)/2;
b = a + delta;

% ----------------------------------------------------------------
% Values of nu to sweep. Include a small value (classic OSWR ~ near
% optimal) up to large values (where optimization should matter a
% lot, per the high-frequency parabolic-like scaling of the KV term)
% ----------------------------------------------------------------
nuList = [0.005, 0.05, 0.5, 1, 2, 5];

nCases = numel(nuList);

% Storage for summary
qOptList   = zeros(nCases,1);
rOptList   = zeros(nCases,1);
rhoOptList = zeros(nCases,1);

errClassicAll = cell(nCases,1);
errOptAll     = cell(nCases,1);

%% ================================================================
% Loop over nu
% ================================================================

for i = 1:nCases

    nu = nuList(i);

    fprintf('==============================================\n');
    fprintf('nu = %.4g\n', nu);
    fprintf('==============================================\n');

    % ------------------------------------------------------------
    % Experiment setup (trapezoid standing wave, Neumann)
    % ------------------------------------------------------------
    uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
    u0_fun = @(x) uExact(x,0);
    v0_fun = @(x) zeros(size(x));
    f_fun  = @(x,t) zeros(size(x));

    % ------------------------------------------------------------
    % Monolithic reference
    % ------------------------------------------------------------
    [~,~,u_mono,~,~,~] = solver( ...
        u0_fun,v0_fun,f_fun,...
        dx,dt,L,T,...
        c,gamma,nu,...
        bcType,...
        1/c,0,...
        @(t)0,@(t)0);

    % ------------------------------------------------------------
    % Classical Robin parameters
    % ------------------------------------------------------------
    qClassic = 1/c;
    rClassic = 0;

    % ------------------------------------------------------------
    % Optimized Robin parameters
    % ------------------------------------------------------------
    [qOpt,rOpt,rhoOpt] = get_optimal_robin_params(c,gamma,nu,a,b,T,dx,dt);

    qOptList(i)   = qOpt;
    rOptList(i)   = rOpt;
    rhoOptList(i) = rhoOpt;

    fprintf('Classical Robin: q = %.6g, r = %.6g\n', qClassic, rClassic);
    fprintf('Optimized Robin: q = %.6g, r = %.6g, predicted rho = %.6e\n\n', ...
        qOpt, rOpt, rhoOpt);

    % ------------------------------------------------------------
    % Classical SWR
    % ------------------------------------------------------------
    [~,~,u_iter_classic,~] = swr_solver( ...
        u0_fun,v0_fun,f_fun,...
        dx,dt,L,T,...
        c,gamma,nu,...
        qClassic,rClassic,...
        delta,iterations);

    % ------------------------------------------------------------
    % Optimized SWR
    % ------------------------------------------------------------
    [~,~,u_iter_opt,~] = swr_solver( ...
        u0_fun,v0_fun,f_fun,...
        dx,dt,L,T,...
        c,gamma,nu,...
        qOpt,rOpt,...
        delta,iterations);

    % ------------------------------------------------------------
    % Errors versus monolithic solution
    % ------------------------------------------------------------
    errClassic = zeros(iterations,1);
    errOpt     = zeros(iterations,1);

    for k = 1:iterations
        errClassic(k) = norm(u_iter_classic(:,end,k)-u_mono(:,end),inf);
        errOpt(k)     = norm(u_iter_opt(:,end,k)-u_mono(:,end),inf);
    end

    errClassicAll{i} = errClassic;
    errOptAll{i}     = errOpt;

end

%% ================================================================
% Plot: one subplot per nu, classic vs optimized convergence
% ================================================================

nCols = ceil(sqrt(nCases));
nRows = ceil(nCases/nCols);

figure('Color','w','Position',[100 100 300*nCols 260*nRows])

for i = 1:nCases

    subplot(nRows,nCols,i)

    semilogy(0:iterations-1, errClassicAll{i}, 'o-', 'LineWidth', 1.5, ...
        'MarkerSize', 3)
    hold on
    semilogy(0:iterations-1, errOptAll{i}, 's-', 'LineWidth', 1.5, ...
        'MarkerSize', 3)

    grid on
    box on

    xlabel('Schwarz iteration')
    ylabel('L_\infty error')

    title(sprintf('\\nu = %.4g', ...
        nuList(i)))

    if i == 1
        legend('Classical (q=1/c, r=0)','Optimized','Location','northeast')
    end

end

sgtitle('SWR convergence: classical vs optimized Robin, across \nu')

%% ================================================================
% Summary plot: predicted rho and optimized (q,r) vs nu
% ================================================================

figure('Color','w','Position',[100 100 900 300])

subplot(1,3,1)
semilogx(nuList, rhoOptList, 'o-', 'LineWidth', 2)
grid on; box on
xlabel('\nu'); ylabel('\rho_{opt}')
title('Predicted optimal convergence factor')

subplot(1,3,2)
semilogx(nuList, qOptList, 'o-', 'LineWidth', 2)
hold on
yline(1/c, '--', 'q_{classic}')
grid on; box on
xlabel('\nu'); ylabel('q^\ast')
title('Optimized q vs \nu')

subplot(1,3,3)
semilogx(nuList, rOptList, 'o-', 'LineWidth', 2)
grid on; box on
xlabel('\nu'); ylabel('r^\ast')
title('Optimized r vs \nu')

sgtitle('Optimized Robin parameters as a function of \nu')

%% ================================================================
% Summary table in the command window
% ================================================================

fprintf('\n=== Summary ===\n');
fprintf('%8s %12s %12s %12s\n','nu','qOpt','rOpt','rhoOpt');
for i = 1:nCases
    fprintf('%8.4g %12.6g %12.6g %12.6e\n', ...
        nuList(i), qOptList(i), rOptList(i), rhoOptList(i));
end

%% ================================================================
% Helper functions (unchanged from single-case script)
% ================================================================

function T = modal_time(k,c,gamma,nu)

    a = gamma + nu*k^2;

    disc = sqrt(a^2 - 4*c^2*k^2);

    lambda1 = (-a + disc)/2;
    lambda2 = (-a - disc)/2;

    A = -lambda2/(lambda1-lambda2);
    B =  lambda1/(lambda1-lambda2);

    T = @(t) A*exp(lambda1*t) + B*exp(lambda2*t);

end

function u = standing_trapezoid_neumann(x,t,M,L,c,gamma,nu)

    u = zeros(size(x));

    for j = 0:M

        m = 2*j + 1;          % odd cosine modes

        k = m*pi/L;

        Tm = modal_time(k,c,gamma,nu);

        u = u + ...
            16/pi^2 * cos(m*pi/4)/m^2 .* ...
            Tm(t) .* cos(k*x);

    end

end