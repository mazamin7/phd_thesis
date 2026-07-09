clear all; close all; clc;
addpath(genpath('../src'));

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
experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005; exp_idx = 4;
% experiment = 'smoothPulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 5;
% experiment = 'trianglePulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0; exp_idx = 6;

% unused
% experiment = 'triangleStandingWaveDirichlet'; T = 2; bcType = 'dirichlet'; gamma = 0.5; nu = 0.005;
% experiment = 'triangleStandingWaveDirichlet'; T = 2; bcType = 'dirichlet'; gamma = 0.5; nu = 0.005;

q = 1/c;
r = 0;

iterations = 25;
delta = L/5;

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
% Optimized Robin parameters
% ------------------------------------------------------------

[qOpt,rOpt,rhoOpt] = get_optimal_robin_params(c,gamma,nu,a,b,T,dx,dt);

fprintf('Optimized Robin:\n');
fprintf('q = %.12g\n',qOpt);
fprintf('r = %.12g\n',rOpt);
fprintf('Predicted convergence factor = %.6e\n\n',rhoOpt);

%% Solve

[x_grid,t_grid,u_iter,v_iter] = swr_solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    qOpt,rOpt,...
    delta,iterations);

u_array = u_iter(:,:,end);
v_array = v_iter(:,:,end);

Nt = length(t_grid);

NxFine = 1000;
dxFine = L/NxFine;

x_fine = ((0:NxFine-1)+0.5)*dxFine;

% Final iterate
u_array_fine = zeros(NxFine,Nt);

for n = 1:Nt
    u_array_fine(:,n) = interp1( ...
        x_grid,...
        u_array(:,n),...
        x_fine,...
        'linear');
end

% Compute monolithic reference solution to calculate SWR iteration error
[~,~,u_mono,~,~,~] = solver( ...
    u0_fun,v0_fun,f_fun,...
    dx,dt,L,T,...
    c,gamma,nu,...
    bcType,q,r,...
    @(t) 0, @(t) 0);

for k = 1:iterations
    err(k) = max(abs(u_iter(:,end,k) - u_mono(:,end)));
end



figure('Color','w')

semilogy(0:iterations-1,err,'o-','LineWidth',2)

grid on
box on

xlabel('Schwarz iteration')
ylabel('L_2 error')
title('Schwarz convergence')

% %% Animate Schwarz iterations
% % All Schwarz iterates for plotting
% u_iter_fine = zeros(NxFine,Nt,iterations);
% 
% for k = 1:iterations
%     for n = 1:Nt
%         u_iter_fine(:,n,k) = interp1( ...
%             x_grid,...
%             u_iter(:,n,k),...
%             x_fine,...
%             'linear');
%     end
% end
% figure('Color','w')
% 
% for k = 1:iterations
% 
%     plot(x_fine,...
%          u_iter_fine(:,end,k),...
%          'b','LineWidth',2)
% 
%     hold on
% 
% 
%     hold off
% 
%     ylim([-1 1])
% 
%     grid on
% 
%     title(sprintf('Iteration %d',k-1))
% 
%     legend('SWR','Exact')
% 
%     drawnow
% 
% end

% %% Animation and Snapshots
% 
% Create a subfolder based on the experiment index
saveFolder = sprintf('snapshots/experiment_%d', exp_idx);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder); 
end
% 
% figure('Color', 'w');
% 
% snapshotTimes = [0.2, 0.5, 1.0];
% 
% % Plot the coarse node positions
% for i = 1:length(x_grid)
%     xline(x_grid(i), 'Color', [1, 0, 0, 0.3], 'LineWidth', 0.5);
% end
% hold on
% 
% % Plot the numerical and exact solutions
% hNum = plot(x_fine, u_array_fine(:,1), 'b', 'LineWidth', 2);
% hExact = plot(x_fine, uExact(x_fine, t_grid(1)), 'r--', 'LineWidth', 2);
% 
% % grid on;
% box on;
% xlabel('x'); ylabel('u');
% ylim([-1 1]);
% hDummy = plot(NaN, NaN, 'Color', [1, 0, 0, 0.5], 'LineWidth', 0.5);
% legend([hNum, hExact, hDummy], 'Numerical', 'Exact', 'Nodes', 'Location', 'best')
% 
% for n = 1:length(t_grid)
%     currentTime = t_grid(n);
% 
%     % Update the plots
%     set(hNum, 'YData', u_array_fine(:,n));
%     set(hExact, 'YData', uExact(x_fine, currentTime));
%     title(sprintf('t = %.4f', currentTime));
%     drawnow;
% 
%     % Check if we need to save a snapshot
%     for s = 1:length(snapshotTimes)
%         if abs(currentTime - snapshotTimes(s)) < (dt/2)
%             filename = fullfile(saveFolder, sprintf('snapshot_t_%.4f.png', currentTime));
%             exportgraphics(gcf, filename, 'Resolution', 300);
%             fprintf('Saved snapshot at t = %.4f\n', currentTime);
%         end
%     end
% end

%% Space-time plot
hSpacetime = figure('Color', 'w');
imagesc(t_grid,x_fine,u_array_fine);
axis xy;
xlabel('Time'); ylabel('Position');
title(sprintf('Space-time solution (Exp %d)', exp_idx));
colorbar;
exportgraphics(hSpacetime, fullfile(saveFolder, 'spacetime_solution.png'), 'Resolution', 300);

%% Surface plot
hSurf = figure('Color', 'w');
surf(t_grid,x_fine,u_array_fine, 'EdgeColor', 'none');
xlabel('Time'); ylabel('Position'); zlabel('u');
title(sprintf('Space-time solution (Exp %d)', exp_idx));
view(45,30); camlight; lighting gouraud;
exportgraphics(hSurf, fullfile(saveFolder, 'surface_solution.png'), 'Resolution', 300);

format long e;
disp('Errors vs monolithic reference:');
disp(err);

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