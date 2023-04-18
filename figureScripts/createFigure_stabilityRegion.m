
% flags
plot_delimiter_lines    = 1;
plot_triangle           = 0;
plot_zero_cover         = 1;
save_figure             = 0;
labels                  = 'on';

% suppress warnings from flags
%#ok<*UNRCH> 

cadence_to_use = 120;
minute_to_second = 1/60;
bpm = cadence_to_use;
metronome_frequency = bpm * 0.5 * minute_to_second; % two beats per metronome cycle
stride_time = metronome_frequency.^(-1);
step_time = stride_time * 0.5;

%% create surf
b_p_min = 0;
b_p_max = 5;
b_v_min = 0;
b_v_max = 1.5;
b_p_number_of_steps = 201; % this value is high for high resolution images, can reduce to speed things up if needed
b_v_number_of_steps = 201;

% set parameters
z_c = 1;
g = 9.81;

% calculate dependent variables
T_total = 1/metronome_frequency * 0.5;
T_step = 1/(2*metronome_frequency);
omega = sqrt(g/z_c);
c = cosh(omega * T_step * 0.5);
s = sinh(omega * T_step * 0.5);

% sample the region of the parameter space where both Eigenvalues are <1 absolute
b_p_grid = linspace(b_p_min, b_p_max, b_p_number_of_steps);
b_v_grid = linspace(b_v_min, b_v_max, b_v_number_of_steps);
steps_per_minute = metronome_frequency * 120;

[B_P,B_V] = meshgrid(b_p_grid,b_v_grid);
Z_rho = zeros(size(B_P));
C_rho = zeros(size(B_P));
alpha_data_rho = zeros(size(B_P));

C_D = zeros(size(B_P));

for i_p = 1 : b_p_number_of_steps
    b_p_here = b_p_grid(i_p);
    for i_v = 1 : b_v_number_of_steps
        b_v_here = b_v_grid(i_v);
        A_here = ...
            [ ...
              (- c*b_p_here + c^2 + s^2),  - (c*b_v_here - 2/omega*c*s); ...
              (2*c - b_p_here)*omega*s, (- b_v_here*omega*s + c^2 + s^2) ...
            ];
        
        % Eigendecomposition to get spectral norm rho(M) (https://en.wikipedia.org/wiki/Convergent_matrix)
        D = eig(A_here);
        rho = max(abs(D));
        Z_rho(i_v, i_p) = rho;
        
        % alternative coloring scheme for poster
        if rho < 1
            C_rho(i_v, i_p) = (1-rho) * 10;
            C_rho(i_v, i_p) = 1-rho;
        else
            C_rho(i_v, i_p) = 0;
            alpha_data_rho(i_v, i_p) = 1;
        end
        
    end
end

% create figure
figure; ax = axes; hold on; title(['step time = ' num2str(T_step) ', Cadence = ' num2str(steps_per_minute) ' steps/min'])
xlabel('b_p'); ylabel('b_v')
xlim([b_p_min, b_p_max]);
ylim([b_v_min, b_v_max]);
if strcmp(labels, 'on')
    file_label = 'withLabels';
elseif strcmp(labels, 'off')
    set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
    set(get(gca, 'xaxis'), 'visible', 'off');
    set(get(gca, 'yaxis'), 'visible', 'off');
    file_label = 'noLabels';
end

% plot surf
surf(B_P, B_V, -Z_rho, 'AlphaData', alpha_data_rho, 'facecolor', 'interp')
clim([-1 0])

shading interp
if strcmp(labels, 'on')
    this_colorbar = colorbar( 'Direction','reverse');
end

% plot delimiter lines
color_triangle = [1 1 1]*0.85;
color_line = [0.5 0.5 0.5];
if plot_delimiter_lines
    % add delimiter lines
    linewidth = 3;
    plot(2*c*[1 1], [b_v_min, b_v_max], 'color', color_line, 'linewidth', linewidth);
    plot([b_p_min, b_p_max], 2*s/omega*[1 1], 'color', color_line, 'linewidth', linewidth);
    
    plot([b_p_min, b_p_max], [b_p_min, b_p_max]*c/(s*omega), 'color', color_line, 'linewidth', linewidth);
end

% plot triangle
if plot_triangle
    P1 = [2*c, 2*s/omega];
    P2 = [2*s^2/c, 2*s/omega];
    P3 = [2*c, 2*c^2/(omega*s)];
    patch([P1(1) P2(1) P3(1)], [P1(2) P2(2) P3(2)], [0 0 0], 'facecolor', color_triangle)
end

% add white plane at zero to cover unstable points
if plot_zero_cover
    xlim([b_p_min, b_p_max]);
    ylim([b_v_min, b_v_max]);
    patch('xdata', [b_p_min b_p_max b_p_max b_p_min], 'ydata', [b_v_min b_v_min b_v_max b_v_max], 'zdata', -1*[1 1 1 1], 'facecolor', 'w', 'edgecolor', 'none')
end

% save
if save_figure
    % save with labels
    filename = ['fig_stabilityRegion_' file_label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
    
    % also export an image that contains the colorbar
    if strcmp(labels, 'off')
        this_colorbar = colorbar('direction', 'reverse', 'ticks', [], 'box', 'off', 'color', 'none');
        print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename '_colorbar'], '-djpeg', '-r600')
        delete(this_colorbar)
    end
end
 






























