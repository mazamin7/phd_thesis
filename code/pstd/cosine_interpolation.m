function u = cosine_interpolation(Um,x,L)

%-------------------------------------------------------------
% interpolate_pstd
%
% Reconstructs the continuous PSTD approximation from the
% nodal values using the DCT-II basis.
%
% Inputs:
%   Um      : DCT coefficients of the function to interpolate
%   x       : evaluation points
%   L       : domain length
%
% Output:
%   u       : interpolated values
%-------------------------------------------------------------

Nold = length(Um);

% Constant mode
u = Um(1)/sqrt(Nold) * ones(size(x));

% Remaining cosine modes
for m = 1:Nold-1
    u = u + ...
        Um(m+1) * ...
        sqrt(2/Nold) * ...
        cos(m*pi*x/L);
end

end