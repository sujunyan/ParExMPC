function plot_sim_results(res_dict)

    % plot the simulation results from different methods stored in res_dict
    % the results are saved in figs/ directory as PNG files

    method_vec = fieldnames(res_dict);
    n_method = length(method_vec);
    plot_options = get_plot_options();

    % Create figs/ directory if it doesn't exist
    if ~exist('figs', 'dir')
        mkdir('figs');
    end

    % Plot state trajectories and save as PNG

    figure('Visible', 'off'); % Create figure without displaying GUI

    lw = 2.0; % line width
    % for method = method_vec
    for i_method = 1:n_method
        method = method_vec{i_method};
        x_vec = res_dict.(method).x_vec;
        nx = size(x_vec, 1);
        for kk = 1:nx
            subplot(nx, 1, kk);
            hold on;
            tspan = 0:(size(x_vec, 2) - 1);
            color = plot_options{i_method}.color;
            ls = plot_options{i_method}.ls;
            plot(tspan, x_vec(kk, :), 'LineWidth', lw, 'DisplayName', method, 'LineStyle',ls, 'Color', color);
            % title(sprintf("State x_%d Trajectory", kk));
            xlabel("Time Step");
            ylabel(sprintf("x_%d", kk));
            grid on;
            legend show;
        end
    end
    hold off;
    saveas(gcf, fullfile('figs', 'state_x_trajectory_admm.png')); % Save as PNG
    close(gcf); % Close the figure

    % Plot running cost and save as PNG
    figure('Visible', 'off'); % Create figure without displaying GUI
    hold on;
    for i_method = 1:n_method
        method = method_vec{i_method};
        J_vec = res_dict.(method).J_vec;
        tspan = 0:(length(J_vec) - 1);
        color = plot_options{i_method}.color;
        ls = plot_options{i_method}.ls;
        plot(tspan, J_vec, 'LineWidth', lw, 'DisplayName', method, 'Color', color, 'LineStyle', ls);
        title("Running Cost J");
        xlabel("Time Step");
        ylabel("Cost");
        grid on;
        legend show;
    end
    hold off;
    saveas(gcf, fullfile('figs', 'running_cost_J_admm.png')); % Save as PNG
    close(gcf); % Close the figure


    % Plot control input and save as PNG
    figure('Visible', 'off'); % Create figure without displaying GUI
    for i_method = 1:n_method
        method = method_vec{i_method};
        u0_vec = res_dict.(method).u0_vec; % Assuming U_vec contains control inputs
        tspan = 0:(size(u0_vec, 2) - 1);
        nu = size(u0_vec, 1);
        for iu = 1:size(u0_vec, 1)
            subplot(nu, 1, iu);
            hold on;
            color = plot_options{i_method}.color;
            ls = plot_options{i_method}.ls;
            plot(tspan, u0_vec(iu, :), 'LineWidth', lw, 'DisplayName', sprintf('%s', method), 'Color', color, 'LineStyle', ls);
            % title("Control Input Trajectory");
            xlabel("Time Step");
            ylabel(sprintf("u^%d_0", iu));
            grid on;
        end

        legend show;
    end
    hold off;
    saveas(gcf, fullfile('figs', 'control_input_trajectory_admm.png')); % Save as PNG
    close(gcf); % Close the figure

end

function plot_options = get_plot_options()
    plot_options = {
        struct('ls', '-', 'color', 'b'), ...
        struct('ls', '--', 'color', 'r'), ...
        struct('ls', ':', 'color', 'g'), ...
        struct('ls', '-.', 'color', 'k') ...
    };
    
end