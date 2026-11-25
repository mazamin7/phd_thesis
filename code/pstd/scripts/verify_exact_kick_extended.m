% Broader checks of the exact kick: zero-viscosity limit and temporal order.
clear; clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));
L=1; c=1; N=32; dx=L/N; T=.1; gamma=.03;
orders=[2 4 6 8]; nus=[0 .05]; steps=[20 40 80];
u0=@(x) exp(-90*(x-.37).^2)+.2*cos(3*pi*x/L);
v0=@(x) .15*sin(2*pi*x/L); f=@(x,t) zeros(size(x));
ns=N/2; x=[((0:ns-1)'+.5)*dx; L/2+((0:ns-1)'+.5)*dx];
lam=-((0:ns-1)'*pi/(L/2)).^2; Q=dct(eye(ns)); LL=idct(diag(lam)*Q); LD=blkdiag(LL,LL);
for nu=nus
  for p=orders
    C=get_residue_matrix(N,p); LT=LD+C/dx^2;
    A=[zeros(N),eye(N); c^2*LT,-gamma*eye(N)+nu*LT];
    yr=expm(T*A)*[u0(x);v0(x)];
    e=zeros(size(steps));
    for k=1:numel(steps)
      dt=T/steps(k);
      [~,~,~,~,u]=ard_solver_exact_kick(u0,v0,f,dx,dt,L,T,c,gamma,nu,'neumannGhost',1,0,p);
      e(k)=norm(u(:,end)-yr(1:N),inf);
    end
    rates=log2(e(1:end-1)./e(2:end));
    fprintf('nu=%g p=%d errors=',nu,p); fprintf(' %.4e',e);
    fprintf(' rates='); fprintf(' %.3f',rates); fprintf('\n');
  end
end
