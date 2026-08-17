for exp_idx = 2:6
    disp(['Running FDTD Exp ', num2str(exp_idx)]);
    cd('code/fdtd/scripts');
    run_solver;
    cd('../../../');
    
    disp(['Running PSTD Exp ', num2str(exp_idx)]);
    cd('code/pstd/scripts');
    run_solver;
    cd('../../../');
end
exit;
