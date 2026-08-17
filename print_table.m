load('table_data.mat');
fprintf('N values: '); disp(N);
fprintf('\n--- FDTD ERRORS ---\n');
for i = 1:size(fdtd_err, 1)
    fprintf('Exp %d: ', i);
    fprintf('%.2e ', fdtd_err(i, :));
    rate = log2(fdtd_err(i, end-1) / fdtd_err(i, end));
    fprintf(' | Rate: %.2f\n', rate);
end

fprintf('\n--- PSTD ERRORS ---\n');
for i = 1:size(pstd_err, 1)
    fprintf('Exp %d: ', i);
    fprintf('%.2e ', pstd_err(i, :));
    rate = log2(pstd_err(i, end-1) / pstd_err(i, end));
    fprintf(' | Rate: %.2f\n', rate);
end
exit;
