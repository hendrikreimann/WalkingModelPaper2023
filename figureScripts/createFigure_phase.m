% simulate SLIP walker with linear step placement feedback

% flags
plot_ap         = 1;
plot_ml         = 1;
save_figures    = 0;
labels          = 'on';

% suppress warnings from flags
%#ok<*UNRCH> 

% set parameters
dt = 1e-4;
metronome_frequency = 1;
T_total = 1;
T_step = 1/(2*metronome_frequency);
z_c = 1;
g = 9.81;
omega = sqrt(g/z_c);
c = cosh(omega * T_step * 0.5);
s = sinh(omega * T_step * 0.5);

% define control parameters - default for cadence 120 (metronome_frequency = 1)
b_o_prg = -0.05;
b_p_prg = 2.5;
b_d_prg = 0.6;

b_o_alt = 0.02;
b_p_alt = 2.5;
b_d_alt = 0.6;

% define initial conditions
p_0_prg = 0.0;
x_0_prg = 0.0;
v_0_prg = b_o_prg*omega / (2*s - b_d_alt*omega);

p_0_alt = 0.0;
x_0_alt = b_o_alt / (2*c - b_p_alt);
v_0_alt = 0;

% define control structs
parameters_prg.control = struct;
parameters_prg.control.T_step = T_step;
parameters_prg.control.b_o = b_o_prg;
parameters_prg.control.b_p = b_p_prg;
parameters_prg.control.b_d = b_d_prg;
parameters_prg.biomechanics = struct;
parameters_prg.biomechanics.z_c = z_c;
parameters_prg.biomechanics.g = g;

parameters_alt.control = struct;
parameters_alt.control.T_step = T_step;
parameters_alt.control.b_o = b_o_alt;
parameters_alt.control.b_p = b_p_alt;
parameters_alt.control.b_d = b_d_alt;
parameters_alt.biomechanics = struct;
parameters_alt.biomechanics.z_c = z_c;
parameters_alt.biomechanics.g = g;

% solve
results_prg = simulateSlipWalker(T_total, parameters_prg, x_0_prg, v_0_prg, p_0_prg, dt, 'progressive');
results_alt = simulateSlipWalker(T_total, parameters_alt, x_0_alt, v_0_alt, p_0_alt, dt, 'alternating');

% find indices for arrows
arrow_interval = 0.05;
arrow_times = 0 : arrow_interval : T_total;
arrows_prg_x = zeros(size(arrow_times));
arrows_prg_v = zeros(size(arrow_times));
arrows_prg_angles = zeros(size(arrow_times));
arrows_prg_dx = zeros(size(arrow_times));
arrows_prg_dv = zeros(size(arrow_times));

arrows_alt_x = zeros(size(arrow_times));
arrows_alt_v = zeros(size(arrow_times));
arrows_alt_angles = zeros(size(arrow_times));
arrows_alt_dx = zeros(size(arrow_times));
arrows_alt_dv = zeros(size(arrow_times));

for i_arrow = 1 : length(arrow_times)
    this_time = arrow_times(i_arrow);
    this_index = findClosestIndex(this_time, results_alt.continuous.time);
    
    arrows_prg_x(i_arrow) = results_prg.continuous.x(this_index);
    arrows_prg_v(i_arrow) = results_prg.continuous.v(this_index);
    arrows_prg_dx(i_arrow) = results_prg.continuous.v(this_index);
    arrows_prg_dv(i_arrow) = results_prg.continuous.a(this_index);

    arrows_alt_x(i_arrow) = results_alt.continuous.x(this_index);
    arrows_alt_v(i_arrow) = results_alt.continuous.v(this_index);
    arrows_alt_dx(i_arrow) = results_alt.continuous.v(this_index);
    arrows_alt_dv(i_arrow) = results_alt.continuous.a(this_index);
end

% separate into parts for steps
step_indices = results_prg.indices.step_indices_continuous;
step_1_indices = 1 : step_indices(1)-1;
step_2_indices = step_indices(1) : step_indices(2)-1;
step_3_indices = step_indices(2) : length(results_prg.continuous.time);
x_step_1_prg = results_prg.continuous.x(step_indices(1));
x_step_2_prg = results_prg.continuous.x(step_indices(2));
x_end_prg = results_prg.continuous.x(end);
p_step_0_prg = results_prg.continuous.p(1);
p_step_1_prg = results_prg.continuous.p(step_indices(1));
p_step_2_prg = results_prg.continuous.p(step_indices(2));

x_step_1_alt = results_alt.continuous.x(step_indices(1));
x_step_2_alt = results_alt.continuous.x(step_indices(2));
x_end_alt = results_alt.continuous.x(end);
p_step_0_alt = results_alt.continuous.p(1);
p_step_1_alt = results_alt.continuous.p(step_indices(1));
p_step_2_alt = results_alt.continuous.p(step_indices(2));

% set visualization parameters
gray = [1 1 1]*0.7;
graphic_specs = struct;
graphic_specs.linewidth = 1;
graphic_specs.orbit_resolution = 100;
graphic_specs.line_color = gray;
graphic_specs.diagonal_color = 'none';
graphic_specs.diagonal_width = 4;
arrow_specs = struct;
arrow_specs.dt = 0.05;
arrow_specs.width = 5;
arrow_specs.length = 8;
arrow_specs.color = gray;

color_orbit = lines(1);
linewidth_orbit = 10;

color_guides = gray;
linewidth_guides = 4;

% ----------------------------------------------------------------------------------------------------------------------
% visualize - anterior-posterior
% ----------------------------------------------------------------------------------------------------------------------

if plot_ap
    figure('position', [0, 0, 800 800], 'renderer', 'painters');
    if strcmp(labels, 'on')
        ax = axes; hold on;
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        ax = axes('position', [0 0 1 1]); hold on
        label = 'noLabels';
    end
    y_min = -1.8;
    y_max = 1.8;
    xlim([0 x_end_prg])
    ylim([y_min y_max])

    x_grid = [0.1 0.2];
    v_grid = [0.5 1 1.5];
%     omega = 3.4715;

    % guidelines
    plot([x_step_1_prg x_step_1_prg], [y_min y_max], 'linewidth', linewidth_guides, 'color', color_guides)
    plot([x_step_2_prg x_step_2_prg], [y_min y_max], 'linewidth', linewidth_guides, 'color', color_guides)

    % orbits
    plot(results_prg.continuous.x(step_1_indices), results_prg.continuous.v(step_1_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);
    plot(results_prg.continuous.x(step_2_indices), results_prg.continuous.v(step_2_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);
    plot(results_prg.continuous.x(step_3_indices), results_prg.continuous.v(step_3_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);    
    
    % phase portraits
    drawSlipPhasePortrait(ax, 0, 0, x_step_1_prg, y_min, y_max, x_grid, [-v_grid v_grid], omega, graphic_specs, arrow_specs)
    drawSlipPhasePortrait(ax, p_step_1_prg, x_step_1_prg, x_step_2_prg, y_min, y_max, p_step_1_prg+[-x_grid x_grid], [-v_grid v_grid], omega, graphic_specs, arrow_specs)
    drawSlipPhasePortrait(ax, p_step_2_prg, x_step_2_prg, x_end_prg, y_min, y_max, p_step_2_prg - x_grid, [-v_grid v_grid], omega, graphic_specs, arrow_specs)
    
    % arrows
    for i_arrow = 1 : length(arrows_prg_x)
        anchor = [arrows_prg_x(i_arrow); arrows_prg_v(i_arrow)]; 
        tip = [arrows_prg_x(i_arrow) + arrows_prg_dx(i_arrow); arrows_prg_v(i_arrow) + arrows_prg_dv(i_arrow)]; 
        plotArrowhead(ax, anchor, tip, arrow_specs.length, arrow_specs.width, arrow_specs.color)
    end

    xlabel('anterior-posterior position', 'fontsize', 14)
    ylabel('anterior-posterior velocity', 'fontsize', 14)
end

if plot_ap && save_figures
    % save with labels
    filename = ['fig_phase_ap_' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')

    close(gcf)
end


% ----------------------------------------------------------------------------------------------------------------------
% visualize - medial-lateral
% ----------------------------------------------------------------------------------------------------------------------

if plot_ml
    fig = figure('position', [0, 0, 800 800], 'renderer', 'painters');
    if strcmp(labels, 'on')
        ax = axes; hold on;
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        ax = axes('position', [0 0 1 1]); hold on
        label = 'noLabels';
    end
    y_min = -0.5;
    y_max = 0.5;
    xlim([0 p_step_1_alt])
    ylim([y_min y_max])

    x_grid = [0.04 0.08 0.12 0.16];
    v_grid = [0.1 0.2 0.3 0.4];
%     omega = 3.4715;

    % guidelines
%     plot([0 p_step_1_alt], [0 0], 'linewidth', linewidth_guides, 'color', color_guides)
    plot([x_step_1_alt x_step_1_alt], [y_min y_max], 'linewidth', linewidth_guides, 'color', color_guides)
%     plot([x_step_2 x_step_2], [y_min y_max], 'linewidth', linewidth_guides, 'color', color_guides)

    % orbits
    plot(results_alt.continuous.x(step_1_indices), results_alt.continuous.v(step_1_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);
    plot(results_alt.continuous.x(step_2_indices), results_alt.continuous.v(step_2_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);
    plot(results_alt.continuous.x(step_3_indices), results_alt.continuous.v(step_3_indices), 'color', color_orbit(1, :), 'linewidth', linewidth_orbit);

    % phase portraits
    drawSlipPhasePortrait(ax, 0, 0, x_step_1_alt, y_min, y_max, x_grid, [-v_grid v_grid], omega, graphic_specs, arrow_specs)
    drawSlipPhasePortrait(ax, p_step_1_alt, x_step_1_alt, p_step_1_alt, y_min, y_max, p_step_1_alt - x_grid, [-v_grid v_grid], omega, graphic_specs, arrow_specs)
%     drawSlipPhasePortrait(ax, p_step_2, x_step_2, x_end, y_min, y_max, p_step_2+[-x_grid], [-v_grid v_grid], omega, graphic_specs, arrow_specs)

    % arrows
    for i_arrow = 1 : length(arrows_alt_x)
        anchor = [arrows_alt_x(i_arrow); arrows_alt_v(i_arrow)]; 
        tip = [arrows_alt_x(i_arrow) + arrows_alt_dx(i_arrow); arrows_alt_v(i_arrow) + arrows_alt_dv(i_arrow)]; 
        plotArrowhead(ax, anchor, tip, arrow_specs.length, arrow_specs.width, arrow_specs.color)
    end

    xlabel('medial-lateral position', 'fontsize', 14)
    ylabel('medial-lateral velocity', 'fontsize', 14)
end

if plot_ml && save_figures
    % save with labels
    filename = ['fig_phase_ml_' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')

    close(gcf)
end

