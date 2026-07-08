clear all; close all; clc;

%% Parameters
c = 1;

% bcType = 'neumannGhost';

for exp_idx = 1:6

    switch exp_idx
        case 1
            experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
        case 2
            experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
        case 3
            experiment = 'standingWave'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005;
        case 4
            experiment = 'trapezoidStandingWaveNeumann'; L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005;
        case 5
            experiment = 'smoothPulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
        case 6
            experiment = 'trianglePulse'; L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
    end

    fprintf('\n=======================================================\n');
    fprintf('Running Experiment %d: %s (Monolithic PSTD)\n', exp_idx, experiment);
    fprintf('=======================================================\n');

    % unused
    q = 1/c;
    r = 0;

    %% Experiment Setup
    switch lower(experiment)
        case 'standingwave'
            k1 = pi/L;
            k3 = 3*pi/L;
            T1 = modal_time(k1,c,gamma,nu);
            T3 = modal_time(k3,c,gamma,nu);
            uExact = @(x,t) cos(k1*x).*T1(t) + 0.3*cos(k3*x).*T3(t);
            u0_fun = @(x) uExact(x,0);
            v0_fun = @(x) zeros(size(x));
            f_fun = @(x,t) zeros(size(x));

        case 'trapezoidstandingwaveneumann'
            uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
            u0_fun = @(x) uExact(x,0);
            v0_fun = @(x) zeros(size(x));
            f_fun  = @(x,t) zeros(size(x));

        case 'smoothpulse'
            if gamma~=0 || nu~=0
                warning(['The manufactured travelling-wave solution is only valid for gamma = nu = 0.']);
            end
            mu    = L/4;
            sigma = L/20;
            M = ceil(c*T/L) + 2;
            uExact = @(x,t) gaussianImages(x,t,mu,c,sigma,L,M);
            u0_fun = @(x) uExact(x,0);
            v0_fun = @(x) zeros(size(x));
            f_fun  = @(x,t) zeros(size(x));

        case 'trianglepulse'
            if gamma~=0 || nu~=0
                warning(['The manufactured travelling-wave solution is only valid for gamma = nu = 0.']);
            end
            mu    = L/4;
            sigma = L/20;
            M = ceil(c*T/L) + 2;
            uExact = @(x,t) triangleImages(x,t,mu,c,sigma,L,M);
            u0_fun = @(x) uExact(x,0);
            v0_fun = @(x) zeros(size(x));
            f_fun  = @(x,t) zeros(size(x));
    end

    %% Refinement study
    N = [25 50 100 200 400 800];

    errProj = zeros(size(N));
    errTrun = zeros(size(N));
    errTot = zeros(size(N));

    for k = 1:length(N)
        dx = L/N(k);
        CFL = 0.9;
        dt = CFL*dx/c;

        [x_grid,t_grid,u_array,v_array,U_array,V_array] = solver( ...
            u0_fun,v0_fun,f_fun,...
            dx,dt,L,T,...
            c,gamma,nu,...
            bcType,q,r);

        Nt = length(t_grid);
        Nx_fine = 10000;
        dx_fine = L/Nx_fine;
        x_fine = ((0:Nx_fine-1)+0.5)*dx_fine;

        u_initial_coarse_num = u_array(:,1);
        u_initial_fine_num = cosine_interpolation(U_array(:,1),x_fine,L);
        u_initial_coarse_ref = uExact(x_grid,t_grid(1));
        u_initial_fine_ref_unfiltered = uExact(x_fine,t_grid(1));

        u_final_coarse_num = u_array(:,end);
        u_final_fine_num = cosine_interpolation(U_array(:,end),x_fine,L);
        u_final_coarse_ref = uExact(x_grid,t_grid(end));
        u_final_fine_ref_unfiltered = uExact(x_fine,t_grid(end));

        U_initial_ref = dct(u_initial_coarse_ref);
        u_initial_fine_ref_filtered = cosine_interpolation(U_initial_ref,x_fine,L);

        U_final_ref = dct(u_final_coarse_ref);
        u_final_fine_ref_filtered = cosine_interpolation(U_final_ref,x_fine,L);
        
        errProj(k) = sqrt(dx_fine)*norm(u_initial_fine_ref_filtered-u_initial_fine_ref_unfiltered);
        errTrun(k) = sqrt(dx)*norm(u_final_coarse_num-u_final_coarse_ref);
        errTot(k) = sqrt(dx_fine)*norm(u_final_fine_num-u_final_fine_ref_unfiltered);
    end

    %% Create output folder
    saveFolder = sprintf('snapshots/experiment_%d', exp_idx);
    if ~exist(saveFolder,'dir')
        mkdir(saveFolder);
    end
    dxVec = L./N;

    %% Discretization Error Plot
    h = figure('Color','w');
    loglog(dxVec,errTrun,'o-','LineWidth',2,'MarkerSize',8,'DisplayName','Measured');
    hold on
    ref = errTrun(1)*(dxVec/dxVec(1)).^2;
    loglog(dxVec,ref,'k--','LineWidth',1.5,'DisplayName','O(\Delta x^2)')
    grid on
    xlabel('\Delta x')
    ylabel('L_2 error')
    title(sprintf('Discretization Error (Monolithic PSTD - Exp %d)', exp_idx))
    legend('Location','best')

    exportgraphics(h, fullfile(saveFolder,'final_nodal_error_monolithic.png'), 'Resolution',300);
    close(h);
end

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