clear all; close all; clc;

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

experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 1;

q = 1/c;
r = 0;

iterations = 100;
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
% Classical Robin parameters
% ------------------------------------------------------------

qClassic = 1/c;
rClassic = 0;

%% ------------------------------------------------------------
% Optimized Robin parameters
% ------------------------------------------------------------

[qOpt,rOpt,rhoOpt] = get_optimal_robin_params(c,gamma,nu,a,b,T,dx,dt);

fprintf('\n');
fprintf('Classical Robin:\n');
fprintf('q = %.12g\n',qClassic);
fprintf('r = %.12g\n\n',rClassic);

fprintf('Optimized Robin:\n');
fprintf('q = %.12g\n',qOpt);
fprintf('r = %.12g\n',rOpt);
fprintf('Predicted convergence factor = %.6e\n\n',rhoOpt);

%% ------------------------------------------------------------
% Classical SWR
% ------------------------------------------------------------

[~,~,u_iter_classic,~] = swr_solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    qClassic,rClassic,...
    delta,iterations);

%% ------------------------------------------------------------
% Optimized SWR
% ------------------------------------------------------------

[~,~,u_iter_opt,~] = swr_solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    qOpt,rOpt,...
    delta,iterations);

%% ------------------------------------------------------------
% Errors versus monolithic solution
% ------------------------------------------------------------

errClassic = zeros(iterations,1);
errOpt     = zeros(iterations,1);

for k = 1:iterations

    errClassic(k) = norm( ...
        u_iter_classic(:,end,k)-u_mono(:,end),inf);

    errOpt(k) = norm( ...
        u_iter_opt(:,end,k)-u_mono(:,end),inf);

end

%% ------------------------------------------------------------
% Plot
% ------------------------------------------------------------

figure('Color','w')

semilogy(0:iterations-1,errClassic,'o-','LineWidth',2)
hold on
semilogy(0:iterations-1,errOpt,'s-','LineWidth',2)

grid on
box on

xlabel('Schwarz iteration')
ylabel('L_\infty error')

legend( ...
    'Robin q=1/c, r=0',...
    'Optimized Robin',...
    'Location','southwest')

title('SWR convergence')

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