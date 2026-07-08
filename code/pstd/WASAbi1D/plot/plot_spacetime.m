function [fig_p, fig_v] = plot_spacetime(t_axis, x_axis, p, v, db_plot, title_str)
%PLOT_SPACETIME Plots the spacetime solutions for pressure and velocity
%   [fig_p, fig_v] = PLOT_SPACETIME(t_axis, x_axis, p, v, title_str) takes
%   in the time axis t_axis, the space axis x_axis, the pressure solution p
%   and the velocity solution v, as well as an optional title string
%   title_str. It then plots the spacetime solutions for both pressure and
%   velocity, with time on the x-axis, space on the y-axis, and pressure or
%   velocity on the z-axis.
%
%   fig_p and fig_v are the figure handles for the pressure and velocity
%   plots, respectively.
%
%   If title_str is provided, it is concatenated before the words
%   'Pressure Solution' or 'Velocity Solution' in the title of the
%   respective plots.

    if db_plot == false
        % Plot p
        fig_p = figure();
        surf(t_axis, x_axis, p,'EdgeColor','none','FaceColor','interp');
        xlabel('Time [s]');
        ylabel('Space [m]');
        zlabel('Pressure');
        if nargin == 6
            title(['Pressure Solution: ', title_str]);
        else
            title('Pressure Solution');
        end
        view(0, 90);  % set view to show from the top
        colorbar;
        
        % Plot v
        fig_v = figure();
        surf(t_axis, x_axis, v,'EdgeColor','none','FaceColor','interp');
        xlabel('Time [s]');
        ylabel('Space [m]');
        zlabel('Velocity');
        if nargin == 6
            title(['Velocity Solution: ', title_str]);
        else
            title('Velocity Solution');
        end
        view(0, 90);  % set view to show from the top
        colorbar;
    else
        % Plot p
        fig_p = figure();
        surf(t_axis, x_axis, max(db(p), -200),'EdgeColor','none','FaceColor','interp');
        hold on;
        z_line = max(db(p(:))) + 10;  % Adjust the value 10 as needed to raise the line
        line([min(t_axis), max(t_axis)], [5, 5], [z_line, z_line], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
        hold off;
        xlabel('Time [s]');
        ylabel('Space [m]');
        zlabel('Pressure (dB)');
        if nargin == 6
            title(['Pressure Solution: ', title_str]);
        else
            title('Pressure Solution');
        end
        view(0, 90);  % set view to show from the top
        colorbar;
        clim([-150 0]);
        
        % Plot v
        fig_v = figure();
        surf(t_axis, x_axis, max(db(v), -200),'EdgeColor','none','FaceColor','interp');
        hold on;
        z_line = max(db(v(:))) + 10;  % Adjust the value 10 as needed to raise the line
        line([min(t_axis), max(t_axis)], [5, 5], [z_line, z_line], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
        hold off;
        xlabel('Time [s]');
        ylabel('Space [m]');
        zlabel('Velocity (dB)');
        if nargin == 6
            title(['Velocity Solution: ', title_str]);
        else
            title('Velocity Solution');
        end
        view(0, 90);  % set view to show from the top
        colorbar;
        clim([-100 10]);
    end
end
