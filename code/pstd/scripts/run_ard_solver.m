clear
close all
clc
addpath(genpath('../src'));
% addpath('WASAbi1D/utils');

%% Parameters
c = 1;

% dx = 0.002;
% dt = 0.002;

% dx = 0.01;
% dt = 0.01;

dx = 0.02;
% dt = 0.02;

% dx = 0.05;
% Use stricter CFL for stability
CFL = 0.2;
dt = CFL * dx;

% dx = 0.1;
% dt = 0.1;

% dx = L/3;
% dt = dx;

% bcType = 'neumannGhost';

% experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 1;
% experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 2;
% experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005; exp_idx = 3;
% experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005; exp_idx = 4;
% experiment = 'smoothPulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 5;
experiment = 'trianglePulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 6;

% unused
q = 1/c;
r = 0;

space_order = 8; % spatial order for interface jump correction (2, 4, 6, or 8)

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
    % Smooth propagating Gaussian pulse
    % ========================================================
    
    case 'smoothpulse'
    
        if gamma~=0 || nu~=0
            warning(['The manufactured travelling-wave solution is only ', ...
                     'valid for gamma = nu = 0.']);
        end
    
        mu    = L/2;
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

%% Solve

[x_grid,x_gridL,x_gridR,t_grid,u,uL,uR,UL,UR] = ard_solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    bcType,q,r,space_order);

Nt = length(t_grid);
Nx = 10000;
dx_fine = L/Nx;
    
if mod(Nx,2)~=0
    error('Nx must be even.');
end
    
NxSub = Nx/2;
Lsub  = L/2;
    
x_fineL = ((0:NxSub-1)' + 0.5)*dx_fine;
x_fineR = Lsub + ((0:NxSub-1)' + 0.5)*dx_fine;
    
x_fine = [x_fineL;x_fineR];


for n = 1:Nt
    u_fine(:,n) = [cosine_interpolation(UL(:,n),x_fineL,Lsub);...
                   cosine_interpolation(UR(:,n),x_fineR - Lsub,Lsub)];
end

%% Animation

figure

hNum = plot(x_fine,u_fine(:,1), ...
    'b','LineWidth',2);

hold on

hExact = plot(x_fine,uExact(x_fine,t_grid(1)), ...
    'r--','LineWidth',2);

grid on
box on

xlabel('x')
ylabel('u')

% ylim([1.1*min(U(:)) 1.1*max(U(:))])
ylim([-1 1])

legend('Numerical','Exact','Location','best')

title(sprintf('t = %.4f',t_grid(1)))

% Compute L2 error against exact solution (for verification)
if strcmpi(experiment,'standingWave')
    uExactVals = uExact(x_fine,t_grid);
    err = sqrt(mean((u_fine - uExactVals).^2,'all'));
    fprintf('L2 error over time: %.3e\n',err);
end

for n = 1:length(t_grid)

    set(hNum,'YData',u_fine(:,n));

    set(hExact,'YData',uExact(x_fine,t_grid(n)));

    title(sprintf('t = %.4f',t_grid(n)));

    drawnow;

end

%% Space-time plot

figure

imagesc(t_grid,x_fine,u_fine)

axis xy

xlabel('Time')
ylabel('Position')

title('Space-time solution')

colorbar

%% Surface plot

figure

surf(t_grid,x_fine,u_fine,...
    'EdgeColor','none')

xlabel('Time')
ylabel('Position')
zlabel('u')

title('Space-time solution')

view(45,30)

camlight
lighting gouraud

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