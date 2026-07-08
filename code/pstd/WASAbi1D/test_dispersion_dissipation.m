clear all, close all, clc;

addpath utils\
addpath sim\
addpath plot\

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;
method = simulation_parameters.method_left;

dh = 4e-4;

% Simulation parameters
if method == 1
    % FDTD 2ord
    lambda_arr = [0.1, 0.2, 0.4, 0.6, 0.8];
    plot_m = 3;
    plot_n = 2;
elseif method == 2
    % FDTD 1ord
    lambda_arr = [0.1, 0.2, 0.4, 0.6, 0.8];
    plot_m = 3;
    plot_n = 2;
elseif method == 3
    % Fourier 2ord
    lambda_arr = [0.5 1 2];
    plot_m = 2;
    plot_n = 2;
elseif method == 4
    % Fourier 1ord
    lambda_arr = [0.5 1 2];
    plot_m = 2;
    plot_n = 2;
end

% Artificial dissipation factor
diss = 1 - eps(1); % - eps(1); % < 1

% Show debug plot?
debug = false;

len_x = test_case_data.len_x;
len_t = test_case_data.len_t;
sigma = len_x/80;

pos_first = len_x/2;
left_first = pos_first - sigma;
right_first = pos_first + sigma;
N = round(right_first/dh) - round(left_first/dh) + 1;

p_gt_fun = test_case_data.p_gt_fun;
x_axis = linspace(left_first, right_first, N);
p_first = p_gt_fun(x_axis, 0)';

f = figure();
position = [100, 100, 500 * plot_n, 200 * plot_m];
set(f, 'Position', position);

subplot(plot_m, plot_n, 1);
hold on;
plot(x_axis, p_first);
xlim([pos_first-sigma, pos_first+sigma]);
ylim([0, 24]);
xlabel("x");
ylabel("p");
title(sprintf("Wave packet at t = %.1f", 0));


fft_size = 1024;
f_max = 1/dh;
f_axis = linspace(-f_max/2,f_max/2-dh,fft_size);

fft_first = fftshift(fft(p_first, fft_size));

f2 = figure();
set(f2, 'Position', position);
subplot(plot_m, plot_n, 1);
plot(f_axis, unwrap(angle(fft_first)));
% xlim([-f_max/2,f_max/2-dh]);
xlim([0, 400]);
% ylim([-360,60]);
xlabel("f");
ylabel(sprintf('\\angle FFT'));
title(sprintf("FFT (unwrapped phase) of the wave packet at t = %.1f", 0));

f3 = figure();
set(f3, 'Position', position);
subplot(plot_m, plot_n, 1);
plot(f_axis, abs(fft_first));
% xlim([-f_max/2,f_max/2-dh]);
xlim([0, 400]);
% ylim([-360,60]);
xlabel("f");
ylabel(sprintf('|FFT|'));
title(sprintf("FFT (magnitude) of the wave packet at t = %.1f", 0));

dct_size = fft_size/2;
dct_axis = linspace(0,f_max/2-dh,dct_size);

f4 = figure();
set(f4, 'Position', position);
subplot(plot_m, plot_n, 1);
plot(dct_axis, dct(p_first, dct_size));
xlim([0, 400]);
% ylim([-360,60]);
xlabel("f");
ylabel(sprintf('DCT'));
title(sprintf("DCT of the wave packet at t = %.1f", 0));


pos_last = len_x/2 + c*1;
left_last = pos_last - sigma;
right_last = pos_last + sigma;

for n = 1:length(lambda_arr)
    lambda = lambda_arr(n);
    dt = dh * lambda / c;
    % Run simulation
    [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, diss, false);

    p_last = p(round(left_last/dh):round(right_last/dh),end);

    str = sprintf("\\lambda = %.2f", lambda);

    figure(f);
    subplot(plot_m, plot_n, 1+n);
    plot(x_axis(round(left_last/dh):round(right_last/dh)), p_last);
    xlim([pos_last-sigma, pos_last+sigma]);
    ylim([0, 12]);
    xlabel("x");
    ylabel("p");
    title(sprintf("Wave packet at t = %.1f - %s", len_t, str));
    
    hold on;
    plot(x_axis(round(left_last/dh):round(right_last/dh)), p_first/2, 'r--');
    legend("Simulated", "Ideal");


    fft_last = fftshift(fft(p_last, fft_size));

    figure(f2);
    subplot(plot_m, plot_n, 1+n);
    plot(f_axis, angle(fft_last./fft_first));
%     xlim([-f_max/2,f_max/2-dh]);
    xlim([0, 400]);
    ylim([-1,1]*pi);
    xlabel("f");
    ylabel(sprintf("\\angle H(f)"));
    yticks([-4 -3 -2 -1 0 1 2 3 4]*pi/4);
    yticklabels({'-4π/4', '-3π/4', '-2π/4', '-π/4', '0', 'π/4', '2π/4', '3π/4', '4π/4'});
    title(sprintf("Phase of the frequency response (FFT ratio) - %s", str));

    hold on;
    plot(f_axis, f_axis*0, 'r--');
    legend("Simulated", "Ideal", "Location", "northwest");

    figure(f3);
    subplot(plot_m, plot_n, 1+n);
    plot(f_axis, 2*abs(fft_last./fft_first));
%     xlim([-f_max/2,f_max/2-dh]);
    xlim([0, 400]);
    ylim([1-0.25,1+0.05]);
    xlabel("f");
    ylabel(sprintf("|H(f)|"));
    title(sprintf("Magnitude of the frequency response (FFT ratio) - %s", str));

    hold on;
    plot(f_axis, f_axis*0 + 1, 'r--');
    legend("Simulated", "Ideal", "Location", "southwest");

    figure(f4);
    subplot(plot_m, plot_n, 1+n);
    plot(dct_axis, dct(p_last, dct_size));
    xlim([0, 400]);
    hold on;
    plot(dct_axis, dct(p_first, dct_size)/2, 'r--');
    xlabel("f");
    ylabel(sprintf("DCT"));
    title(sprintf("DCT of the wave packet at t = %.1f - %s", len_t, str));
    legend("Simulated", "Ideal");
end








