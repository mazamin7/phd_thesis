t_grid = 0:0.1:1;
vals = rand(11, 10); % 11 times, 10 spaces
v = interp1(t_grid, vals, 0.5);
disp(size(v))
