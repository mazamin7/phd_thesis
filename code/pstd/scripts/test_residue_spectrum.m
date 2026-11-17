% TEST_RESIDUE_SPECTRUM Verify spectral properties of the ARD interface residual.
% The Kick uses (nu/dx^2)*C_res, where C_res is returned by
% get_residue_matrix. Since nu/dx^2 is positive, signs of its eigenvalues
% are the same as those of C_res.

clear; clc;
addpath(genpath('../src'));

orders = [2, 4, 6, 8];
Nxs = [20, 40, 80, 160, 320, 640, 800];
tolImag = 1e-10;
tolPos = 1e-10;

rows = [];
for Nx = Nxs
    for order = orders
        C_res = full(get_residue_matrix(Nx, order));
        lambda = eig(C_res);
        scale = max(1, max(abs(lambda)));
        maxImag = max(abs(imag(lambda)));
        maxReal = max(real(lambda));
        minReal = min(real(lambda));
        realWithinTol = maxImag <= tolImag * scale;
        nonPositiveWithinTol = maxReal <= tolPos * scale;
        rows(end+1,:) = [Nx, order, maxImag, minReal, maxReal, ...
            realWithinTol, nonPositiveWithinTol]; %#ok<SAGROW>
        fprintf('Nx=%4d order=%d max|Im(lambda)|=%.3e minRe=%.6e maxRe=%.6e real=%d nonpositive=%d\n', ...
            Nx, order, maxImag, minReal, maxReal, realWithinTol, nonPositiveWithinTol);
    end
end

results = array2table(rows, 'VariableNames', ...
    {'Nx','Order','MaxAbsImag','MinReal','MaxReal','RealWithinTol','NonPositiveWithinTol'});
disp(results);
writetable(results, 'residue_spectrum_results.csv');

if all(results.RealWithinTol)
    fprintf('PASS: all tested C_res matrices have real eigenvalues within the stated tolerance.\n');
else
    fprintf('FAIL: some tested C_res matrices have non-real eigenvalues beyond tolerance.\n');
end
if all(results.NonPositiveWithinTol)
    fprintf('PASS: all tested C_res matrices are non-positive within the stated tolerance.\n');
else
    fprintf('FAIL: some tested C_res matrices have positive eigenvalues beyond tolerance.\n');
end
