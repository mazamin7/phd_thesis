function plot_snapshot(x_axis,len_x,p,v,c0,f,db_plot)
%PLOT_SNAPSHOT Summary of this function goes here
%   Detailed explanation goes here
    small_font = true;

    % Plot
    figure(f);

    if db_plot == false
        % Plot p
        subplot(2,1,1);
        plot(x_axis, p);
        title('Pressure');
        xlim([0,len_x]);
        ylim([-1,1]*2e-1);
        xlabel("x");
        ylabel("p");
    
        % Plot v
        subplot(2,1,2);
        plot(x_axis, v);
        title('Velocity');
        xlim([0,len_x]);
        ylim([-c0,c0]*5e-1);
        xlabel("x");
        ylabel("v");
    else
        % Plot p
        subplot(211);
        plot(x_axis, db(p));
        hold on;
        line([5 5], [-150 0], 'Color', 'red', 'LineStyle', '--');
        hold off;
        title('Pressure (dB)');
        grid on;
        xlim([0,len_x]);
        ylim([-140 10]);
        if small_font == true
            yticks(-140:10:10);
        else
            yticks(-140:30:10);

            ax = gca;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.Title.FontSize = 20;
        end
        xlabel("x");
        ylabel("p (dB)");
    
        % Plot v
        subplot(212);
        plot(x_axis, db(v));
        hold on;
        line([5 5], [-150 0], 'Color', 'red', 'LineStyle', '--');
        hold off;
        title('Velocity (dB)');
        grid on;
        xlim([0,len_x]);
        ylim([-130 20]);
        if small_font == true
            yticks(-140:10:10);
        else
            yticks(-140:30:10);

            ax = gca;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.Title.FontSize = 20;
        end
        xlabel("x");
        ylabel("v (dB)");
    end
end

