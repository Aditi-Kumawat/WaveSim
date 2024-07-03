% Define radius vector
r_vect = 0.1:0.1:10;

%-------------------------------------------------------%%
soil_medium = fns_inptMatPara.select_soil_medium();
disp(['selectedRfFldr: ', soil_medium])
%-------------------------------------------------------%%
d_J = 2.5; % Depth value
%-------------------------------------------------------%%
dir = fns_inptMatPara.form_dir_Iparas(soil_medium, d_J);
% Load the displacement data
load(fullfile(dir, 'Uphi_t_mat.mat'));
load(fullfile(dir, 'Ur_t_mat.mat'));
load(fullfile(dir, 'Uz_t_mat.mat'));
Uphi_t_mat = Uphi_t_mat * 1e3; % Convert to mm
Ur_t_mat = Ur_t_mat * 1e3; % Convert to mm
Uz_t_mat = Uz_t_mat * 1e3; % Convert to mm

% Load the velocity data
load(fullfile(dir, 'Vphi_t_mat.mat'));
load(fullfile(dir, 'Vr_t_mat.mat'));
load(fullfile(dir, 'Vz_t_mat.mat'));
Vphi_t_mat = Vphi_t_mat * 1e3; % Convert to mm/s
Vr_t_mat = Vr_t_mat * 1e3; % Convert to mm/s
Vz_t_mat = Vz_t_mat * 1e3; % Convert to mm/s

% Define the time vector
t_vect = 0:0.01:6;

% Define colors
ha_cl=@colors;
cl = {ha_cl('ball blue'), ha_cl('carmine red'), ha_cl('black')};

% Create a figure for the animation
figure('Position', [300, 300, 1000, 200]); % Width x Height

% Displacement Animation
% animate_component(r_vect, t_vect, Ur_t_mat, Uphi_t_mat, Uz_t_mat, cl, ...
%     'Displacement $(\mathrm{mm})$', {'$u_r$', '$u_\phi$', '$u_z$'}, 'Displacement components with time', 'displacement_animation.mp4');

% Velocity Animation
animate_component(r_vect, t_vect, Vr_t_mat, Vphi_t_mat, Vz_t_mat, cl, ...
    'Velocity $(\mathrm{mm/s})$', {'$v_r$', '$v_\phi$', '$v_z$'}, 'Velocity Components with Time', 'velocity_animation.mp4');

function animate_component(r_vect, t_vect, Ur_t_mat, Uphi_t_mat, Uz_t_mat, cl, y_label, legend_labels, title_text, filename)
    % Calculate the magnitude of the displacement or velocity vector
    Ur_max = max(max(real(Ur_t_mat)));
    Ur_min = min(min(real(Ur_t_mat)));
    Uphi_max = max(max(real(Uphi_t_mat)));
    Uphi_min = min(min(real(Uphi_t_mat)));
    Uz_max = max(max(real(Uz_t_mat)));
    Uz_min = min(min(real(Uz_t_mat)));

    % Create the initial plot with legends
    h1 = plot(r_vect, real(Ur_t_mat(:, 1)), 'Color', cl{1}, 'LineWidth', 2); hold on;
    h2 = plot(r_vect, real(Uphi_t_mat(:, 1)), 'Color', cl{2}, 'LineWidth', 2);
    h3 = plot(r_vect, real(Uz_t_mat(:, 1)), 'Color', cl{3}, 'LineWidth', 1.5);
    hold off;

    % Set axis limits
    xlim([min(r_vect) max(r_vect)]);
    ylim([min([Ur_min, Uphi_min, Uz_min]), max([Ur_max, Uphi_max, Uz_max])]);

    % Add labels, title, and legend using LaTeX
    xlabel('$r\, (\mathrm{km})$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(y_label, 'Interpreter', 'latex', 'FontSize', 12);
    legend(legend_labels, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 14,'Box','off');
    title(title_text, 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    set(gca, 'Color', 'white');
    % Create a VideoWriter object
    v = VideoWriter(filename, 'MPEG-4');
    v.FrameRate = 10; % Set frame rate
    open(v);

    % Loop through time steps to create the animation
    for i = 1:length(t_vect)
        % Update the data for the plots
        set(h1, 'YData', real(Ur_t_mat(:, i)));
        set(h2, 'YData', real(Uphi_t_mat(:, i)));
        set(h3, 'YData', real(Uz_t_mat(:, i)));

        % Update the title with the current time using sprintf for consistent formatting
        title([title_text, ', $t = ', sprintf('%.2f', t_vect(i)), '\, \mathrm{s}$'], 'Interpreter', 'latex');

        % Capture the plot as a frame
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    % Close the video file
    close(v);
end
