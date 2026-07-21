function [x_grid,x_gridL,x_gridR,t_grid,u,uL,uR,UL_array,UR_array] = ard_solver(u0_fun,v0_fun,f_fun,dx,dt,L,T,c,gamma,nu,bcType,q,r,space_order)
    
    if nargin < 14
        space_order = 2;
    end
    Nx = round(L/dx);
    
    if mod(Nx,2)~=0
        error('Nx must be even.');
    end
    
    NxSub = Nx/2;
    Lsub  = L/2;
    
    % -------------------------------------------------------------
    % Non-overlapped sub-domains (standard ARD)
    % -------------------------------------------------------------
    x_gridL = ((0:NxSub-1)'+0.5)*dx;               % left sub-domain nodes
    x_gridR = Lsub + ((0:NxSub-1)'+0.5)*dx;       % right sub-domain nodes

    % Original grids for output/storage (same as physical grids)
    x_gridL_no = x_gridL;
    x_gridR_no = x_gridR;

    x_grid = [x_gridL_no; x_gridR_no];
    
    t_grid = 0:dt:T;
    
    Nt = length(t_grid);
    
    %-------------------------------------------------------------
    % Storage
    %-------------------------------------------------------------
    
    uL = zeros(NxSub,Nt);
    uR = zeros(NxSub,Nt);
    UL_array = zeros(NxSub,Nt);
    UR_array = zeros(NxSub,Nt);
    
    u  = zeros(Nx,Nt);
    
    %-------------------------------------------------------------
    % Initial conditions
    %-------------------------------------------------------------
    
    % Initial conditions on overlapped grids
    u0L_full = u0_fun(x_gridL);
    u0R_full = u0_fun(x_gridR);
    
    % Store only non-overlapped portion for output arrays
    uL(:,1) = u0L_full;
    uR(:,1) = u0R_full;
    
    % Combine into global solution (non-overlapped ordering)
    u(:,1) = [uL(:,1); uR(:,1)];
    
    % Initial modal coefficients from the full overlapped fields
    UL = dct(u0L_full);
    VL = dct(v0_fun(x_gridL));
    UR = dct(u0R_full);
    VR = dct(v0_fun(x_gridR));
    %-------------------------------------------------------------
    % No interface forcing needed for standard ARD – continuity will be enforced after modal propagation.
    
    % (Placeholder – nothing to compute here)
    
    % Compute modal eigenvalues (negative squared wave numbers)
    m = (0:NxSub-1)';
    k = m*pi/Lsub;
    lambda = -k.^2;
    
    [S,Tm] = modal_propagators( ...
        lambda,...
        dt,...
        c,...
        gamma,...
        nu);
    
    %-------------------------------------------------------------
    % Penalty method parameters
    %-------------------------------------------------------------
    
    % No explicit phi basis needed; interface correction will be applied directly to modal vectors.
        
    %-------------------------------------------------------------
    % Initial interface forces
    %-------------------------------------------------------------
    
    C_res = get_residue_matrix(Nx, space_order);
    
    [FcorrL, FcorrR, uL_phys, uR_phys] = get_interface_forces(UL, UR, VL, VR, NxSub, C_res, dx, c, nu);

    %-------------------------------------------------------------
    % Time stepping
    %-------------------------------------------------------------
    
    print_step = max(1, floor(Nt / 10));
    
    for n = 1:Nt-1
        if mod(n, print_step) == 0
            fprintf('  ard_solver progress: %3.0f%%\n', (n/Nt)*100);
        end
        
        t = t_grid(n);
        
        % Physical forcing (zero in most tests)
        fL = f_fun(x_gridL,t);
        fR = f_fun(x_gridR,t);
        FmL = dct(fL);
        FmR = dct(fR);
        
        % First half-kick with interface forces
        VL = VL + (dt/2) * FcorrL;
        VR = VR + (dt/2) * FcorrR;
        
        % Drift: modal propagation with physical forcing
        [UL,VL] = modal_step(UL, VL, FmL, S, Tm);
        [UR,VR] = modal_step(UR, VR, FmR, S, Tm);
        
        % Recompute interface forces
        [FcorrL, FcorrR, uL_phys, uR_phys] = get_interface_forces(UL, UR, VL, VR, NxSub, C_res, dx, c, nu);
        
        % Second half-kick
        VL = VL + (dt/2) * FcorrL;
        VR = VR + (dt/2) * FcorrR;
        
        % Store results
        uL(:,n+1) = uL_phys;
        uR(:,n+1) = uR_phys;
        
        UL_array(:,n+1) = UL;
        UR_array(:,n+1) = UR;
        
        % Concatenate to global solution
        u(:,n+1) = [uL(:,n+1); uR(:,n+1)];
    end
    
    % End of function
    
    return

end


function [Um,Vm] = modal_step( ...
    Um,Vm,...
    Fm,...
    S,T)

    Nm = length(Um);
    
    for j = 1:Nm
    
        z = [Um(j);Vm(j)];
    
        z = S{j}*z + T{j}*[0;Fm(j)];
    
        Um(j) = z(1);
        Vm(j) = z(2);
    
    end
end

function [FcorrL, FcorrR, uL_phys, uR_phys] = get_interface_forces(UL, UR, VL, VR, NxSub, C_res, dx, c, nu)
        uL_phys = idct(UL);
        uR_phys = idct(UR);
        vL_phys = idct(VL);
        vR_phys = idct(VR);
        
        u_phys = [uL_phys; uR_phys];
        v_phys = [vL_phys; vR_phys];
        fCorr_global = (c^2/dx^2) * (C_res * u_phys) + (nu/dx^2) * (C_res * v_phys);
        fCorrL = fCorr_global(1:NxSub);
        fCorrR = fCorr_global(NxSub+1:end);
        
        FcorrL = dct(fCorrL);
        FcorrR = dct(fCorrR);
    end