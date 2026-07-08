clear all, close all, clc;

addpath utils\
addpath sim\
addpath plot\

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Plotting ground truth
plot_ground_truth(test_case_data, false);


% Simulation parameters
CFL = 0.8;
dh = 1e-2;
dt = dh * CFL / c;

% Artificial dissipation factor
diss = 1 - eps(1);

% Show debug plot?
debug = true;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, diss, false);

% Plotting simulation
[fig_p, fig_v] = plot_spacetime(t_axis,x_axis,p,v,false,'Simulation');

% Save simulation as figures and animation
save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v);






