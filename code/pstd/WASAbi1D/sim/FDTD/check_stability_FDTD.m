function [stable] = check_stability_FDTD(c, dt, dh, space_order, diss, asymptotic)
%
% Only works for space_order == 6
%

    assert(space_order == 2 || space_order == 4 || space_order == 6 || space_order == 8);
    n = space_order/2;

    if n == 1
        coefs = [1, -2];
    elseif n == 2
        coefs = [-1/12, 4/3, -5/2];
    elseif n == 3
        coefs = [1/90, -3/20, 3/2, -49/18];
    else
        coefs = [-1/560, 8/315, -1/5, 8/5, -205/72];
    end

    var = 1/2 * (sum(abs(coefs(1:end-1))) - coefs(end)/2);

    CFL = 1 / sqrt(var);

    if asymptotic
        stable = dt < dh / c * CFL; % CFL condition
    else
        stable = dt <= dh / c * CFL;
    end

    stable = stable && diss < 1;

end