function [stable] = check_stability(c, dt, dh, space_order, fourier, diss)

    % check stability
    if fourier == false % FDTD
        stable = check_stability_FDTD(c, dt, dh, space_order, diss, true);
    elseif fourier == true % Fourier
        cond1 = check_stability_Fourier(diss);
        cond2 = check_stability_FDTD(c, dt, dh, space_order, diss, false);
        stable = cond1 && cond2;
    end
end