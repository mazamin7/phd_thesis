clear all; close all; clc;
addpath(genpath('../src'));

%% Parameters
c = 1;

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
experiment = 'trianglePulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 6;

% unused
% experiment = 'triangleStandingWaveDirichlet'; L = 1; T = 2; bcType = 'dirichlet'; gamma = 0; nu = 0; exp_idx = 7;
% experiment = 'triangleStandingWaveDirichlet'; L = 1; T = 2; bcType = 'dirichlet'; gamma = 0.5; nu = 0.005; exp_idx = 8;

q = 1/c;
r = 0;

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
    
        mu    = L/2;
        sigma = L/20;
    
        M = ceil(c*T/L) + 2;
    
        uExact = @(x,t) triangleImages(x,t,mu,c,sigma,L,M);
    
        u0_fun = @(x) uExact(x,0);
        v0_fun = @(x) zeros(size(x));
        f_fun  = @(x,t) zeros(size(x));

    otherwise

        error('Unknown experiment.');

end

%% Refinement study

% N = [10 25 50 100 200 400];
N = [25 50 100 200 400 800];
% N = [100 200 400 800 1600];

% err = zeros(size(N));

errProj = zeros(size(N));
errTrun = zeros(size(N));
errTot = zeros(size(N));

for k = 1:length(N)

    dx = L/N(k);

    CFL = 1;

    dt = CFL*dx/c;

    [x_grid,t_grid,u_array,v_array] = solver( ...
        u0_fun,v0_fun,f_fun,...
        dx,dt,L,T,...
        c,gamma,nu,...
        bcType,q,r);

    Nt = length(t_grid);
    Nx = 10000;
    
    dx_fine = L/Nx;
    
    x_fine = linspace(0,L,Nx+1).';
    
    % u_array_fine = interp1( ...
    %     x_grid, ...
    %     u_array(:,end), ...
    %     x_fine, ...
    %     'linear');
    % 
    % u_ref = uExact(x_fine,t_grid(end));
    
    % err(k) = sqrt(dx_fine)*norm(u_array_fine-u_ref);

    u_initial_coarse_num = u_array(:,1);
    u_initial_fine_num = interp1(x_grid,u_initial_coarse_num,x_fine,'linear');
    u_initial_coarse_ref = uExact(x_grid,t_grid(1));
    u_initial_fine_ref_unfiltered = uExact(x_fine,t_grid(1));

    u_final_coarse_num = u_array(:,end);
    u_final_fine_num = interp1(x_grid,u_final_coarse_num,x_fine,'linear');
    u_final_coarse_ref = uExact(x_grid,t_grid(end));
    u_final_fine_ref_unfiltered = uExact(x_fine,t_grid(end));

    % err(k) = sqrt(dx_fine)*norm(u_array_fine-u_ref);

    % Filter the initial condition
    u_initial_fine_ref_filtered = interp1(x_grid,u_initial_coarse_ref,x_fine,'linear');

    % Filter the exact solution
    u_final_fine_ref_filtered = interp1(x_grid,u_final_coarse_ref,x_fine,'linear');
    
    errProj(k) = sqrt(dx_fine)*norm(u_initial_fine_ref_filtered-u_initial_fine_ref_unfiltered);
    errTrun(k) = sqrt(dx)*norm(u_final_coarse_num-u_final_coarse_ref);
    errTot(k) = sqrt(dx_fine)*norm(u_final_fine_num-u_final_fine_ref_unfiltered);

end

% %% Experimental order
% 
% eoc = log(err(1:end-1)./err(2:end))./log(2);
% 
% disp(table( ...
%     N(:),...
%     err(:),...
%     [NaN;eoc(:)],...
%     'VariableNames',{'N','Error','EOC'}))
% 
% %% Convergence Plot
% % Create a subfolder based on the experiment index
% saveFolder = sprintf('snapshots/experiment_%d', exp_idx);
% if ~exist(saveFolder, 'dir')
%     mkdir(saveFolder); 
% end
% 
% hConv = figure('Color', 'w');
% % Set markers and colors
% loglog(L./N, err, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical Error');
% hold on;
% % Reference line for O(dx^2)
% loglog(L./N, (L./N).^2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(\Delta x^2)');
% 
% grid on;
% set(gca, 'FontSize', 12, 'LineWidth', 1.2); % Professional axis settings
% xlabel('Spatial Step \Delta x', 'FontSize', 14);
% ylabel('L_2 Error', 'FontSize', 14);
% title('Convergence Analysis', 'FontSize', 16);
% legend('Location', 'best', 'FontSize', 12);
% 
% % Save to the same folder as your previous snapshots
% exportgraphics(hConv, fullfile(saveFolder, 'convergence_study.png'), 'Resolution', 300);

% %% Experimental orders
% 
% eocProj = log(errProj(1:end-1)./errProj(2:end))/log(2);
% eocTrun = log(errTrun(1:end-1)./errTrun(2:end))/log(2);
% eocTot  = log(errTot (1:end-1)./errTot (2:end))/log(2);
% 
% disp(table( ...
%     N(:), ...
%     errProj(:), ...
%     [NaN; eocProj(:)], ...
%     errTrun(:), ...
%     [NaN; eocTrun(:)], ...
%     errTot(:), ...
%     [NaN; eocTot(:)], ...
%     'VariableNames', ...
%     {'N','Initial Projection Error','Initial Projection Error OC', ...
%      'Final Nodal Error','Final Nodal Error OC', ...
%      'Final Continuous Error','Final Continuous Error OC'}))

%% Create output folder

saveFolder = sprintf('snapshots/experiment_%d', exp_idx);

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

dxVec = L./N;

% %% --------------------------------------------------------
% % Initial Projection Error
% % --------------------------------------------------------
% 
% h = figure('Color','w');
% 
% loglog(dxVec,errProj,'o-','LineWidth',2,...
%     'MarkerSize',8,...
%     'DisplayName','Measured');
% hold on
% 
% ref = errProj(1)*(dxVec/dxVec(1)).^2;
% 
% loglog(dxVec,ref,'k--','LineWidth',1.5,...
%     'DisplayName','O(\Delta x^2)')
% 
% grid on
% xlabel('\Delta x')
% ylabel('L_2 error')
% title('Initial Projection Error')
% legend('Location','best')
% 
% exportgraphics(h,...
%     fullfile(saveFolder,'initial_projection_error.png'),...
%     'Resolution',300);

%% --------------------------------------------------------
% Final Nodal Error
% --------------------------------------------------------

h = figure('Color','w');

loglog(dxVec,errTrun,'o-','LineWidth',2,...
    'MarkerSize',8,...
    'DisplayName','Measured');
hold on

ref = errTrun(1)*(dxVec/dxVec(1)).^2;

loglog(dxVec,ref,'k--','LineWidth',1.5,...
    'DisplayName','O(\Delta x^2)')

grid on
xlabel('\Delta x')
ylabel('L_2 error')
title('Final Nodal Error')
legend('Location','best')

exportgraphics(h,...
    fullfile(saveFolder,'final_nodal_error.png'),...
    'Resolution',300);

% %% --------------------------------------------------------
% % Final Continuous Error
% % --------------------------------------------------------
% 
% h = figure('Color','w');
% 
% loglog(dxVec,errTot,'o-','LineWidth',2,...
%     'MarkerSize',8,...
%     'DisplayName','Measured');
% hold on
% 
% ref = errTot(1)*(dxVec/dxVec(1)).^2;
% 
% loglog(dxVec,ref,'k--','LineWidth',1.5,...
%     'DisplayName','O(\Delta x^2)')
% 
% grid on
% xlabel('\Delta x')
% ylabel('L_2 error')
% title('Final Continuous Error')
% legend('Location','best')
% 
% exportgraphics(h,...
%     fullfile(saveFolder,'final_continuous_error.png'),...
%     'Resolution',300);

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