addpath('code/fdtd/scripts');
addpath('code/fdtd/src');
addpath('code/pstd/scripts');
addpath('code/pstd/src');

N = [25 50 100 200 400 800];
fdtd_err = zeros(6, 6);
pstd_err = zeros(6, 6);

c = 1;
for exp_idx = 1:6
    disp(['Running Exp ', num2str(exp_idx)]);
    switch exp_idx
        case 1
            L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
            k1 = pi/L; k3 = 3*pi/L;
            T1 = modal_time(k1,c,gamma,nu); T3 = modal_time(k3,c,gamma,nu);
            uExact = @(x,t) cos(k1*x).*T1(t) + 0.3*cos(k3*x).*T3(t);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun = @(x,t) zeros(size(x));
        case 2
            L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
            uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun  = @(x,t) zeros(size(x));
        case 3
            L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005;
            k1 = pi/L; k3 = 3*pi/L;
            T1 = modal_time(k1,c,gamma,nu); T3 = modal_time(k3,c,gamma,nu);
            uExact = @(x,t) cos(k1*x).*T1(t) + 0.3*cos(k3*x).*T3(t);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun = @(x,t) zeros(size(x));
        case 4
            L = 1; T = 2; bcType = 'neumannGhost'; gamma = 0.5; nu = 0.005;
            uExact = @(x,t) standing_trapezoid_neumann(x,t,1000,L,c,gamma,nu);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun  = @(x,t) zeros(size(x));
        case 5
            L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
            mu = L/4; sigma = L/20; M = ceil(c*T/L) + 2;
            uExact = @(x,t) gaussianImages(x,t,mu,c,sigma,L,M);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun  = @(x,t) zeros(size(x));
        case 6
            L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
            mu = L/4; sigma = L/20; M = ceil(c*T/L) + 2;
            uExact = @(x,t) triangleImages(x,t,mu,c,sigma,L,M);
            u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun  = @(x,t) zeros(size(x));
    end
    
    q = 1/c; r = 0;
    
    for k = 1:length(N)
        dx = L/N(k);
        % FDTD
        dt_f = 0.9 * dx / c;
        addpath('code/fdtd/src');
        [~,~,u_array_f,~,~,~] = solver(u0_fun,v0_fun,f_fun,dx,dt_f,L,T,c,gamma,nu,bcType,q,r);
        rmpath('code/fdtd/src');
        
        x_nodes_f = linspace(0, L, N(k)+1)';
        exact_vals_f = reshape(uExact(x_nodes_f, T), [], 1);
        err_f = u_array_f(:,end) - exact_vals_f;
        fdtd_err(exp_idx, k) = sqrt(dx * sum(err_f.^2));
        
        % PSTD
        dt_p = 0.9*dx/c;
        addpath('code/pstd/src');
        [~,~,u_array_p,~,~,~] = solver(u0_fun,v0_fun,f_fun,dx,dt_p,L,T,c,gamma,nu,bcType,q,r);
        rmpath('code/pstd/src');
        
        x_nodes_p = ((0:N(k)-1)+0.5)'*dx;
        exact_vals_p = reshape(uExact(x_nodes_p, T), [], 1);
        err_p = u_array_p(:,end) - exact_vals_p;
        pstd_err(exp_idx, k) = sqrt(dx * sum(err_p.^2));
    end
end
save('table_data.mat', 'fdtd_err', 'pstd_err', 'N');
disp('Data saved.');
exit;

function T = modal_time(k,c,gamma,nu)
    a = gamma + nu*k.^2;
    disc = sqrt(complex(a.^2 - 4*c^2*k.^2));
    lambda1 = (-a + disc)/2;
    lambda2 = (-a - disc)/2;
    A = -lambda2./(lambda1-lambda2);
    B =  lambda1./(lambda1-lambda2);
    T = @(t) A.*exp(lambda1*t) + B.*exp(lambda2*t);
end

function u = standing_trapezoid_neumann(x,t,M,L,c,gamma,nu)
    u = zeros(size(x));
    for j = 0:M
        m = 2*j + 1;
        k = m*pi/L;
        Tm = modal_time(k,c,gamma,nu);
        u = u + 16/pi^2 * cos(m*pi/4)/m^2 .* Tm(t) .* cos(k*x);
    end
end

function u = gaussianImages(x,t,mu,c,sigma,L,M)
    u = zeros(size(x));
    for m = -M:M
        u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)-c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)+c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)-c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)+c*t)/sigma).^2);
    end
end

function u = triangleImages(x,t,mu,c,sigma,L,M)
    triangle = @(z) max(1-abs(z)/sigma,0);
    u = zeros(size(x));
    for m = -M:M
        u = u + 0.5*triangle(x-(2*m*L + mu)-c*t);
        u = u + 0.5*triangle(x-(2*m*L + mu)+c*t);
        u = u + 0.5*triangle(x-(2*m*L - mu)-c*t);
        u = u + 0.5*triangle(x-(2*m*L - mu)+c*t);
    end
end
