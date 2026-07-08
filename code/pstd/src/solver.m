function [x_grid,t_grid,u_array,v_array,U_array,V_array] = solver( ...
    u0_fun,v0_fun,f_fun,dx,dt,L,T,...
    c,gamma,nu,...
    bcType,q,r)

%-------------------------------------------------------------
% Cell-centered grid
%-------------------------------------------------------------

Nx = round(L/dx);

x_grid = ((0:Nx-1)' + 0.5)*dx;

t_grid = 0:dt:T;

Nt = length(t_grid);

u_array = zeros(Nx,Nt);
v_array = zeros(Nx,Nt);
U_array = zeros(Nx,Nt);
V_array = zeros(Nx,Nt);

%-------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------

u0 = u0_fun(x_grid);
v0 = v0_fun(x_grid);

u_array(:,1) = u0;
v_array(:,1) = v0;

%-------------------------------------------------------------
% Modal transform
%-------------------------------------------------------------

Um = dct(u0);
Vm = dct(v0);

U_array(:,1) = Um;
V_array(:,1) = Vm;

%-------------------------------------------------------------
% Continuous eigenvalues
%-------------------------------------------------------------

m = (0:Nx-1)';

k = m*pi/L;

lambda = -k.^2;

%-------------------------------------------------------------
% Time stepping
%-------------------------------------------------------------

for n = 1:Nt-1

    t = t_grid(n);

    Fm = dct(f_fun(x_grid,t));

    for j = 1:Nx

        d = gamma - nu*lambda(j);

        %-------------------------------------------------
        % Four damping regimes
        %-------------------------------------------------

        if lambda(j)==0 && d==0

            S = [1 dt;
                 0 1];

            Tm = [dt dt^2/2;
                  0 dt];

        elseif lambda(j)==0

            S = [ ...
                1 (1-exp(-d*dt))/d;
                0 exp(-d*dt)];

            Tm = [ ...
                dt (exp(-d*dt)-1+d*dt)/d^2;
                0 (1-exp(-d*dt))/d];

        else

            discr = d^2 + 4*c^2*lambda(j);

            if abs(discr) < 1e-12

                e = exp(-d*dt/2);

                S = e * [ ...
                    1+d*dt/2 dt;
                    -d^2*dt/4 1-d*dt/2];

            elseif discr < 0

                omega = sqrt(-discr)/2;

                e = exp(-d*dt/2);

                S = e*( ...
                    cos(omega*dt)*eye(2) + ...
                    sin(omega*dt)*[ ...
                    d/(2*omega) 1/omega;
                    c^2*lambda(j)/omega -d/(2*omega)]);

            else

                delta = sqrt(discr);

                e = exp(-d*dt/2);

                S = e*( ...
                    cosh(delta*dt/2)*eye(2) + ...
                    sinh(delta*dt/2)*[ ...
                    d/delta 2/delta;
                    2*c^2*lambda(j)/delta -d/delta]);

            end

            A = [ ...
                0 1;
                c^2*lambda(j) -d];

            Tm = A\(S-eye(2));

        end

        z = [Um(j);Vm(j)];

        z = S*z + Tm*[0;Fm(j)];

        Um(j) = z(1);
        Vm(j) = z(2);

    end

    u_array(:,n+1) = idct(Um);
    v_array(:,n+1) = idct(Vm);
    U_array(:,n+1) = Um;
    V_array(:,n+1) = Vm;

end

end