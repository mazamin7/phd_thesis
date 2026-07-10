addpath(genpath('../src'));
% clear all; close all; clc;

%% Parameters
c = 1;

% dx = 0.001;
% dt = 0.001;

% dx = 0.01;
% dt = 0.01;

dx = 0.05;
dt = 0.05;

% bcType = 'dirichlet';
% bcType = 'neumannGhost';
% bcType = 'robinClassic';
% bcType = 'robinBdf';
% bcType = 'robinGhost';

% experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 1;
% experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 2;
% experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005; exp_idx = 3;
% experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005; exp_idx = 4;
% experiment = 'smoothPulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 5;
% experiment = 'trianglePulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 6;

% unused
% experiment = 'triangleStandingWaveDirichlet'; T = 2; bcType = 'dirichlet'; gamma = 0.5; nu = 0.005;
% experiment = 'triangleStandingWaveDirichlet'; T = 2; bcType = 'dirichlet'; gamma = 0.5; nu = 0.005;

experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 2;

q = 1/c;
r = 0;

iterations = 30;
iterationsToMeasure = 30;
delta = L/10;

a = (L-delta)/2;
b = a + delta;

%% Experiment

switch lower(experiment)

    % ========================================================
    % Standing wave (exactly represented by the cosine basis)
    % ========================================================

    case 'standingwave'

        k1 = pi/L;
        k3 = 3*pi/L;

        T1 = modal_time(k1,c,gamma,nu);
        T3 = modal_time(k3,c,gamma,nu);

        uExact = @(x,t) ...
            cos(k1*x).*T1(t) ...
          + 0.3*cos(k3*x).*T3(t);

        u0_fun = @(x) uExact(x,0);

        v0_fun = @(x) zeros(size(x));

        f_fun = @(x,t) zeros(size(x));

    % ========================================================
    % Non-smooth standing wave (Neumann)
    % ========================================================

    case 'trapezoidstandingwaveneumann'

        uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
        
        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun  = @(x,t) zeros(size(x));

    % ========================================================
    % Non-smooth standing wave (Dirichlet)
    % ========================================================

    case 'trianglestandingwavedirichlet'

        uExact = @(x,t) standing_triangle_dirichlet(x,t,1000,L,c,gamma,nu);

        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun = @(x,t) zeros(size(x));

    % ========================================================
    % Smooth propagating Gaussian pulse
    % ========================================================
    
    case 'smoothpulse'
    
        if gamma~=0 || nu~=0
            warning(['The manufactured travelling-wave solution is only ', ...
                     'valid for gamma = nu = 0.']);
        end
    
        mu    = L/4;
        sigma = L/20;
    
        % number of image cells (increase if T is larger)
        M = ceil(c*T/L) + 2;
    
        uExact = @(x,t) gaussianImages(x,t,mu,c,sigma,L,M);
    
        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun  = @(x,t) zeros(size(x));

    % ========================================================
    % Non-smooth propagating triangular pulse
    % ========================================================
    
    case 'trianglepulse'
    
        if gamma~=0 || nu~=0
            warning(['The manufactured travelling-wave solution is only ', ...
                     'valid for gamma = nu = 0.']);
        end
    
        mu    = L/4;
        sigma = L/20;
    
        M = ceil(c*T/L) + 2;
    
        uExact = @(x,t) triangleImages(x,t,mu,c,sigma,L,M);
    
        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun  = @(x,t) zeros(size(x));

    otherwise

        error('Unknown experiment.');

end

%% ------------------------------------------------------------
% Monolithic reference
% ------------------------------------------------------------

[~,~,u_mono,~,~,~] = solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    bcType,...
    1/c,0,...
    @(t)0,@(t)0);

%% ------------------------------------------------------------
% Parameter grid
% ------------------------------------------------------------

qVals = linspace(0,1.5/c,21);
rVals = linspace(0,2,21);

Err = zeros(length(rVals),length(qVals));
Rho = zeros(length(rVals),length(qVals));

wmin = pi/T;
wmax = pi/dt;
omega = logspace(log10(wmin),log10(wmax),4000);

%% ------------------------------------------------------------
% Sweep
% ------------------------------------------------------------

for iq = 1:length(qVals)

    fprintf('%d / %d\n',iq,length(qVals));

    for ir = 1:length(rVals)

        q = qVals(iq);
        r = rVals(ir);

        %-----------------------------------------
        % Continuous convergence factor
        %-----------------------------------------

        k = sqrt((omega.^2-1i*gamma*omega) ...
                ./(c^2+1i*nu*omega));

        Lambda = r + 1i*q*omega;

        rho = abs((1i*k.*sinh(1i*k*a)-Lambda.*cosh(1i*k*a))./(1i*k.*sinh(1i*k*b)+Lambda.*cosh(1i*k*b))).^2 ...
             .*exp(-2*imag(k)*delta);

        Rho(ir,iq) = max(rho);

        %-----------------------------------------
        % SWR
        %-----------------------------------------

        [~,~,u_iter,~] = swr_solver( ...
            u0_fun,v0_fun,f_fun,...
            dx,dt,L,T,...
            c,gamma,nu,...
            q,r,...
            delta,iterationsToMeasure);

        Err(ir,iq) = norm( ...
            u_iter(:,end,iterationsToMeasure)-u_mono(:,end),inf);

    end

end

%% ------------------------------------------------------------
% Locate minima
% ------------------------------------------------------------

[minErr,idx] = min(Err(:));
[irErr,iqErr] = ind2sub(size(Err),idx);

[minRho,idx] = min(Rho(:));
[irRho,iqRho] = ind2sub(size(Rho),idx);

fprintf('\n');
fprintf('Measured optimum\n');
fprintf('q = %.6f\n',qVals(iqErr));
fprintf('r = %.6f\n',rVals(irErr));

fprintf('\n');
fprintf('Continuous optimum\n');
fprintf('q = %.6f\n',qVals(iqRho));
fprintf('r = %.6f\n',rVals(irRho));

%% ------------------------------------------------------------
% Measured error
% ------------------------------------------------------------

fig1 = figure('Color','w');

imagesc(qVals,rVals,log10(Err))

axis xy
colorbar

xlabel('q')
ylabel('r')

title(sprintf('log_{10} error after %d Schwarz iterations', ...
    iterationsToMeasure))

hold on

plot(qVals(iqErr),rVals(irErr),...
    'wo',...
    'MarkerSize',12,...
    'LineWidth',2)

plot(qVals(iqRho),rVals(irRho),...
    'r+',...
    'MarkerSize',14,...
    'LineWidth',2)

legend('Measured optimum','Continuous optimum')

if ~exist('snapshots/swr_parameters', 'dir')
    mkdir('snapshots/swr_parameters');
end
exportgraphics(fig1, 'snapshots/swr_parameters/swr_empirical_error.png', 'Resolution', 300);

%% ------------------------------------------------------------
% Predicted convergence factor
% ------------------------------------------------------------

fig2 = figure('Color','w');

imagesc(qVals,rVals,Rho)

axis xy
colorbar

xlabel('q')
ylabel('r')

title('Predicted max convergence factor')

hold on

plot(qVals(iqErr),rVals(irErr),...
    'wo',...
    'MarkerSize',12,...
    'LineWidth',2)

plot(qVals(iqRho),rVals(irRho),...
    'r+',...
    'MarkerSize',14,...
    'LineWidth',2)

legend('Measured optimum','Continuous optimum')

exportgraphics(fig2, 'snapshots/swr_parameters/swr_predicted_rho.png', 'Resolution', 300);

function T = modal_time(k,c,gamma,nu)

    a = gamma + nu*k^2;
    
    disc = sqrt(a^2 - 4*c^2*k^2);
    
    lambda1 = (-a + disc)/2;
    lambda2 = (-a - disc)/2;
    
    A = -lambda2/(lambda1-lambda2);
    B =  lambda1/(lambda1-lambda2);
    
    T = @(t) A*exp(lambda1*t) + B*exp(lambda2*t);

end

function u = standing_triangle_dirichlet(x,t,M,L,c,gamma,nu)

u = zeros(size(x));

for j = 0:M

    m = 2*j + 1;

    k = m*pi/L;

    Tm = modal_time(k,c,gamma,nu);

    u = u + ...
        8/pi^2 * (-1)^j/m^2 * ...
        Tm(t) .* sin(k*x);

end

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

function u = gaussianImages(x,t,mu,c,sigma,L,M)

u = zeros(size(x));

for m = -M:M

    % image at +mu
    u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)-c*t)/sigma).^2);
    u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)+c*t)/sigma).^2);

    % image at -mu
    u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)-c*t)/sigma).^2);
    u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)+c*t)/sigma).^2);

end

end

function u = triangleImages(x,t,mu,c,sigma,L,M)

triangle = @(z) max(1-abs(z)/sigma,0);

u = zeros(size(x));

for m = -M:M

    % image at +mu
    u = u + 0.5*triangle(x-(2*m*L + mu)-c*t);
    u = u + 0.5*triangle(x-(2*m*L + mu)+c*t);

    % image at -mu
    u = u + 0.5*triangle(x-(2*m*L - mu)-c*t);
    u = u + 0.5*triangle(x-(2*m*L - mu)+c*t);

end

end