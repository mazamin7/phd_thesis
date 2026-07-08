function [S,T] = modal_propagators( ...
    lambda,dt,c,gamma,nu)

    Nm = length(lambda);
    
    S = cell(Nm,1);
    T = cell(Nm,1);
    
    for j = 1:Nm
    
        d = gamma - nu*lambda(j);
    
        if lambda(j)==0 && d==0
    
            S{j} = [1 dt;
                    0 1];
    
            T{j} = [dt dt^2/2;
                    0 dt];
    
        elseif lambda(j)==0
    
            S{j} = ...
                [1 (1-exp(-d*dt))/d;
                 0 exp(-d*dt)];
    
            T{j} = ...
                [dt (exp(-d*dt)-1+d*dt)/d^2;
                 0 (1-exp(-d*dt))/d];
    
        else
    
            discr = d^2 + 4*c^2*lambda(j);
    
            if abs(discr)<1e-12
    
                e = exp(-d*dt/2);
    
                S{j} = e * ...
                    [1+d*dt/2 dt;
                     -d^2*dt/4 1-d*dt/2];
    
            elseif discr<0
    
                omega = sqrt(-discr)/2;
    
                e = exp(-d*dt/2);
    
                S{j} = e*( ...
                    cos(omega*dt)*eye(2) + ...
                    sin(omega*dt)* ...
                    [ d/(2*omega)          1/omega;
                      c^2*lambda(j)/omega -d/(2*omega)]);
    
            else
    
                delta = sqrt(discr);
    
                e = exp(-d*dt/2);
    
                S{j} = e*( ...
                    cosh(delta*dt/2)*eye(2) + ...
                    sinh(delta*dt/2)* ...
                    [ d/delta             2/delta;
                      2*c^2*lambda(j)/delta -d/delta]);
    
            end
    
            A = [ ...
                0 1;
                c^2*lambda(j) -d];
    
            T{j} = A\(S{j}-eye(2));
    
        end
    
    end

end
