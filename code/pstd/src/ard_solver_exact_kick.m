function [x_grid,x_gridL,x_gridR,t_grid,u,uL,uR,UL_array,UR_array] = ard_solver_exact_kick( ...
    u0_fun,v0_fun,f_fun,dx,dt,L,T,c,gamma,nu,bcType,q,r,space_order)
%ARD_SOLVER_EXACT_KICK ARD with an exact interface-kick subflow.
%   The interface residual acts only on a small set of rows. Its exponential
%   is therefore applied through a small matrix function, without forming a
%   dense exponential of the full-grid matrix.

    if nargin < 14 || isempty(space_order)
        space_order = 2;
    end

    Nx = round(L/dx);
    if mod(Nx,2) ~= 0
        error('Nx must be even.');
    end

    NxSub = Nx/2;
    Lsub = L/2;
    x_gridL = ((0:NxSub-1)'+0.5)*dx;
    x_gridR = Lsub + ((0:NxSub-1)'+0.5)*dx;
    x_grid = [x_gridL; x_gridR];
    t_grid = 0:dt:T;
    Nt = length(t_grid);

    uL = zeros(NxSub,Nt);
    uR = zeros(NxSub,Nt);
    u = zeros(Nx,Nt);
    UL_array = zeros(NxSub,Nt);
    UR_array = zeros(NxSub,Nt);

    uL(:,1) = u0_fun(x_gridL);
    uR(:,1) = u0_fun(x_gridR);
    u(:,1) = [uL(:,1); uR(:,1)];

    UL = dct(uL(:,1));
    VL = dct(v0_fun(x_gridL));
    UR = dct(uR(:,1));
    VR = dct(v0_fun(x_gridR));
    UL_array(:,1) = UL;
    UR_array(:,1) = UR;

    modes = (0:NxSub-1)';
    lambda = -(modes*pi/Lsub).^2;
    [S,Tm] = modal_propagators(lambda,dt,c,gamma,nu);

    % Factor C_res = U*V, where U inserts values at affected rows and
    % V = C_res(rows,:). Then V*U = C_res(rows,rows), whose dimension is at
    % most space_order for the current 1D residual stencil.
    C_res = get_residue_matrix(Nx,space_order);
    rows = find(any(C_res,2));
    C_rows = C_res(rows,:);
    M = full(C_res(rows,rows));

    tau = dt/2;
    kick_a = nu/dx^2;
    kick_b = c^2/dx^2;
    z = tau*kick_a;
    Phi1 = matrix_phi1(z*M);

    print_step = max(1,floor(Nt/10));
    for n = 1:Nt-1
        if mod(n,print_step) == 0
            fprintf('  ard_solver_exact_kick progress: %3.0f%%\n',(n/Nt)*100);
        end

        t = t_grid(n);

        % First exact half-kick.
        [VL,VR] = apply_exact_kick( ...
            UL,UR,VL,VR,NxSub,rows,C_rows,Phi1,z,tau*kick_b);

        % Exact modal drift. Midpoint sampling of f preserves second order
        % for smooth time-dependent forcing.
        fL = f_fun(x_gridL,t+dt/2);
        fR = f_fun(x_gridR,t+dt/2);
        [UL,VL] = modal_step(UL,VL,dct(fL),S,Tm);
        [UR,VR] = modal_step(UR,VR,dct(fR),S,Tm);

        % Second exact half-kick.
        [VL,VR] = apply_exact_kick( ...
            UL,UR,VL,VR,NxSub,rows,C_rows,Phi1,z,tau*kick_b);

        uL(:,n+1) = idct(UL);
        uR(:,n+1) = idct(UR);
        u(:,n+1) = [uL(:,n+1);uR(:,n+1)];
        UL_array(:,n+1) = UL;
        UR_array(:,n+1) = UR;
    end
end

function [VL_new,VR_new] = apply_exact_kick( ...
    UL,UR,VL,VR,NxSub,rows,C_rows,Phi1,z,forcing_scale)
    u_phys = [idct(UL);idct(UR)];
    v_phys = [idct(VL);idct(VR)];

    % Exact kick for v' = a*C_res*v + b*C_res*u with u fixed. For
    % C_res=U*V, M=V*U, z=tau*a, and Phi1(X)=(exp(X)-I)/X, the update is
    %   v_new = v + U*(z*Phi1(z*M)*V*v + tau*b*Phi1(z*M)*V*u).
    % Phi1 is evaluated without dividing by X, including the nu=0 limit.
    delta_v = z*Phi1*(C_rows*v_phys) + ...
              forcing_scale*Phi1*(C_rows*u_phys);
    v_phys(rows) = v_phys(rows) + delta_v;

    VL_new = dct(v_phys(1:NxSub));
    VR_new = dct(v_phys(NxSub+1:end));
end

function P = matrix_phi1(X)
% Evaluate phi_1(X) from one small block exponential, robust at X=0.
    r = size(X,1);
    if r == 0
        P = zeros(0,0);
        return
    end
    X = full(X);
    B = [X, eye(r); zeros(r), zeros(r)];
    E = expm(B);
    P = E(1:r,r+1:2*r);
end

function [Um,Vm] = modal_step(Um,Vm,Fm,S,T)
    for j = 1:length(Um)
        z = S{j}*[Um(j);Vm(j)] + T{j}*[0;Fm(j)];
        Um(j) = z(1);
        Vm(j) = z(2);
    end
end
