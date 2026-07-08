function [qOpt, rOpt, rhoOpt] = get_optimal_robin_params(c, gamma, nu, a, b, T, dx, dt)
% GET_OPTIMAL_ROBIN_PARAMS  Optimize Robin transmission parameters (q,r)
% for OSWR applied to the damped wave equation (Neumann boundary case).
%
% Inputs:
%   c, gamma, nu : physical parameters (wave speed, friction, viscosity)
%   a, b         : interface locations of the overlapping subdomains
%   T            : simulation time window
%   dx, dt       : mesh size and time step
%
% Outputs:
%   qOpt, rOpt   : optimized Robin parameters
%   rhoOpt       : optimal (min-max) convergence factor

    % ------------------------------------------------------------
    % Frequency range: bounded by simulation window and time step
    % ------------------------------------------------------------
    wMin  = pi / T;
    wMax  = pi / dt;
    omega = logspace(log10(wMin), log10(wMax), 4000);

    % ------------------------------------------------------------
    % Initial guess and bounds
    % ------------------------------------------------------------
    x0 = [1/c, gamma/(2*c)];   % [q0, r0]
    lb = [0, 0];               % positivity constraint

    % ------------------------------------------------------------
    % Optimization settings
    % ------------------------------------------------------------
    options = optimoptions('patternsearch', ...
        'Display',          'iter', ...
        'UseCompletePoll',  true, ...
        'UseCompleteSearch',true, ...
        'MeshTolerance',    1e-10);

    objectiveFun = @(x) contConvergenceObjective(x(1), x(2), omega, c, gamma, nu, a, b);
    % objectiveFun = @(x) discConvergenceObjective(x(1), x(2), omega, c, gamma, nu, a, b, dx, dt);

    [x, rhoOpt] = patternsearch(objectiveFun, x0, [], [], [], [], lb, [], options);

    qOpt = x(1);
    rOpt = x(2);
end

% ==================================================================
function value = contConvergenceObjective(q, r, omega, c, gamma, nu, a, b)
% Continuous (semi-discrete in space/time) convergence factor, Neumann case.

    k      = sqrt((omega.^2 - 1i*gamma*omega) ./ (c^2 + 1i*nu*omega));
    Lambda = r + 1i*q*omega;

    num = Lambda .* cosh(1i*k*a) - 1i*k .* sinh(1i*k*a);
    den = Lambda .* cosh(1i*k*b) + 1i*k .* sinh(1i*k*b);

    rho   = abs(num ./ den).^2;
    value = max(rho);
end

% ==================================================================
function value = discConvergenceObjective(q, r, omega, c, gamma, nu, a, b, dx, dt)
% Fully discrete convergence factor, Neumann case.

    A   = cos(omega*dt) + 1i*(gamma*dt/2) .* sin(omega*dt);
    den = 2*(c*dt/dx)^2 + 1i*2*nu*dt/dx^2 .* sin(omega*dt);
    s   = (1 - A) ./ den;

    k = (2/dx) * asin(sqrt(s));
    k(imag(k) > 0) = -k(imag(k) > 0);   % enforce correct branch

    Lambda = r + 1i*q*omega;

    num = Lambda .* cosh(1i*k*a) - 1i*k .* sinh(1i*k*a);
    den = Lambda .* cosh(1i*k*b) + 1i*k .* sinh(1i*k*b);

    rho   = abs(num ./ den).^2;
    value = max(rho);
end