addpath('code/fdtd/scripts');
addpath('code/fdtd/src');
L = 2; T = 2; bcType = 'neumannGhost'; gamma = 0; nu = 0;
c = 1; mu = L/4; sigma = L/20; M = ceil(c*T/L) + 2;
uExact = @(x,t) gaussianImages(x,t,mu,c,sigma,L,M);
u0_fun = @(x) uExact(x,0); v0_fun = @(x) zeros(size(x)); f_fun  = @(x,t) zeros(size(x));
q = 1; r = 0; dx = L/25; dt_f = dx/c;
cd code/fdtd/src
[~,~,u_array_f,~,~,~] = solver(u0_fun,v0_fun,f_fun,dx,dt_f,L,T,c,gamma,nu,bcType,q,r);
cd ../../../
disp(['u_array_f size: ', num2str(size(u_array_f(:,end)))]);

function u = gaussianImages(x,t,mu,c,sigma,L,M)
    u = zeros(size(x));
    for m = -M:M
        u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)-c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L + mu)+c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)-c*t)/sigma).^2);
        u = u + 0.5*exp(-0.5*((x-( 2*m*L - mu)+c*t)/sigma).^2);
    end
end
