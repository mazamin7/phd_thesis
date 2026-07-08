function generate_sweeps()
% Function to run the SWR parameter sweeps and save the output figures

close all; clc;

% Set up output directory
out_dir = '../../latex/Images/fdtd_snapshots/swr_parameters/';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% ------------------------------------------------------------
% Run Part 2 (Viscosity Sweep)
% ------------------------------------------------------------
fprintf('Running convergence_test_swr_part2...\n');
% Because convergence_test_swr_part2 has "clear all", it will clear variables
% but since we are in a function, it only clears the function workspace? No, 
% `clear all` in a script called from a function behaves weirdly. 
% Better to just redefine out_dir after running it, or use absolute paths.
run('convergence_test_swr_part2.m');

figs = findall(0, 'Type', 'figure');
[~, sortIdx] = sort([figs.Number]);
figs = figs(sortIdx);

if length(figs) >= 2
    saveas(figs(1), '../../latex/Images/fdtd_snapshots/swr_parameters/swr_nu_sweep.png');
    saveas(figs(2), '../../latex/Images/fdtd_snapshots/swr_parameters/swr_optimized_params.png');
end
close all;

% ------------------------------------------------------------
% Run Part 3 (Empirical vs Theoretical Optimum)
% ------------------------------------------------------------
fprintf('Running convergence_test_swr_part3...\n');
run('convergence_test_swr_part3.m');

figs = findall(0, 'Type', 'figure');
[~, sortIdx] = sort([figs.Number]);
figs = figs(sortIdx);

if length(figs) >= 2
    saveas(figs(1), '../../latex/Images/fdtd_snapshots/swr_parameters/swr_empirical_error.png');
    saveas(figs(2), '../../latex/Images/fdtd_snapshots/swr_parameters/swr_predicted_rho.png');
end
close all;

fprintf('All figures generated successfully.\n');
end
