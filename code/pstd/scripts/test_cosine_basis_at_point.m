clear
close all
clc
addpath(genpath('../src'));

%% Parameters

L = 1;

N = 32;

dx = L/N;

x = ((0:N-1)' + 0.5)*dx;

%% Compare analytical basis with MATLAB DCT

err = zeros(N,1);

for i = 1:N

    %----------------------------------------------
    % Analytical modal coefficients of delta(x-x_i)
    %----------------------------------------------

    phi = cosine_basis_at_point(x(i),L,N);

    %----------------------------------------------
    % Discrete delta
    %----------------------------------------------

    f = zeros(N,1);
    f(i) = 1/dx;

    %----------------------------------------------
    % Modal coefficients from MATLAB
    %----------------------------------------------

    F = dx*dct(f);

    %----------------------------------------------
    % Error
    %----------------------------------------------

    err(i) = norm(phi-F,inf);

end

fprintf('Maximum error = %.3e\n',max(err));

figure
semilogy(1:N,err,'o-','LineWidth',2)
grid on
xlabel('Grid point')
ylabel('Infinity norm error')
title('Consistency of cosine\_basis\_at\_point')