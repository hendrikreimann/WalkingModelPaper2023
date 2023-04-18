% Create graphs for Figure 1: Example Trajectories

% flags
individual_figures          = 0;
save_figures                = 0;
labels                      = 'off';

% suppress warnings from flags
%#ok<*UNRCH> 

% visualization parameters
figure_width_traj = 600;
figure_height_traj = 300;
figure_width_path = 600;
figure_height_path = 950;

% colors
color_r = hex2rgb('#009ADE');
color_1 = hex2rgb('#FFC61E');
color_2 = hex2rgb('#DE400B');
color_3 = hex2rgb('#920BDB');

% set parameters
dt = 1e-4;
metronome_frequency = 1;
T_total = 3;
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
p_0_prg_1 = 0.0;
x_0_prg_1 = 0.0;
v_0_prg_1 = b_o_prg*omega / (2*s - b_d_alt*omega);

p_0_alt_1 = 0.0;
x_0_alt_1 = b_o_alt / (2*c - b_p_alt);
v_0_alt_1 = 0;

p_0_prg_2 = p_0_prg_1 - 0.1;
x_0_prg_2 = x_0_prg_1;
v_0_prg_2 = v_0_prg_1 + 0.1;

p_0_alt_2 = p_0_alt_1 - 0.1;
x_0_alt_2 = x_0_alt_1;
v_0_alt_2 = v_0_alt_1 + 0.1;

p_0_prg_3 = p_0_prg_1 + 0.1;
x_0_prg_3 = x_0_prg_1;
v_0_prg_3 = v_0_prg_1 - 0.1;

p_0_alt_3 = p_0_alt_1 + 0.1;
x_0_alt_3 = x_0_alt_1;
v_0_alt_3 = v_0_alt_1 - 0.1;

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
results_prg_1 = simulateSlipWalker(T_total, parameters_prg, x_0_prg_1, v_0_prg_1, p_0_prg_1, dt, 'progressive');
results_alt_1 = simulateSlipWalker(T_total, parameters_alt, x_0_alt_1, v_0_alt_1, p_0_alt_1, dt, 'alternating');
results_prg_2 = simulateSlipWalker(T_total, parameters_prg, x_0_prg_2, v_0_prg_2, p_0_prg_2, dt, 'progressive');
results_alt_2 = simulateSlipWalker(T_total, parameters_alt, x_0_alt_2, v_0_alt_2, p_0_alt_2, dt, 'alternating');
results_prg_3 = simulateSlipWalker(T_total, parameters_prg, x_0_prg_3, v_0_prg_3, p_0_prg_3, dt, 'progressive');
results_alt_3 = simulateSlipWalker(T_total, parameters_alt, x_0_alt_3, v_0_alt_3, p_0_alt_3, dt, 'alternating');

% visualize
colors = lines(3);
markersize_cop = 12;
linewidth_cop = 2;
linewidth_com = 8;

if ~individual_figures
    fig = figure('position', [0, 0, 1600 800]);
    tiledlayout(3, 3, 'TileSpacing', 'Compact', 'padding', 'compact');
end

% anterior-posterior, position
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(1); hold on
end
plot(results_prg_1.continuous.time, results_prg_1.continuous.x, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_1.continuous.time, results_prg_1.continuous.p, 'color', color_r, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
plot(results_prg_2.continuous.time, results_prg_2.continuous.x, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_2.continuous.time, results_prg_2.continuous.p, 'color', color_1, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
plot(results_prg_3.continuous.time, results_prg_3.continuous.x, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_3.continuous.time, results_prg_3.continuous.p, 'color', color_2, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
xlabel('time (s)', 'fontsize', 14)
ylabel('anterior-posterior position', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ap_pos' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% anterior-posterior, velocity
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(4); hold on
end
plot(results_prg_1.continuous.time, results_prg_1.continuous.v, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_2.continuous.time, results_prg_2.continuous.v, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_3.continuous.time, results_prg_3.continuous.v, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
xlabel('time (s)', 'fontsize', 14)
ylabel('anterior-posterior velocity', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ap_vel' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% anterior-posterior, relative position
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(7); hold on
end
plot(results_prg_1.continuous.time, results_prg_1.continuous.q, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_2.continuous.time, results_prg_2.continuous.q, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_prg_3.continuous.time, results_prg_3.continuous.q, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
xlabel('time (s)', 'fontsize', 14)
ylabel('anterior-posterior relative position', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ap_posrel' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% medial-lateral, position
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(2); hold on
end
plot(results_alt_1.continuous.time, results_alt_1.continuous.x, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_1.continuous.time, results_alt_1.continuous.p, 'color', color_r, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
plot(results_alt_2.continuous.time, results_alt_2.continuous.x, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_2.continuous.time, results_alt_2.continuous.p, 'color', color_1, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
plot(results_alt_3.continuous.time, results_alt_3.continuous.x, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_3.continuous.time, results_alt_3.continuous.p, 'color', color_2, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral position', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ml_pos' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% medial-lateral, velocity
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(5); hold on
end
plot(results_alt_1.continuous.time, results_alt_1.continuous.v, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_2.continuous.time, results_alt_2.continuous.v, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_3.continuous.time, results_alt_3.continuous.v, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral velocity', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ml_vel' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% medial-lateral, relative position
if individual_figures
    figure('position', [100 100 figure_width_traj figure_height_traj]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(8); hold on
end
plot(results_alt_1.continuous.time, results_alt_1.continuous.q, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_2.continuous.time, results_alt_2.continuous.q, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_3.continuous.time, results_alt_3.continuous.q, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral relative position', 'fontsize', 14)
if individual_figures && save_figures
    filename = ['fig_example_ml_posrel' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end

% path
if individual_figures
    figure('position', [100 100 figure_width_path figure_height_path]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(3, [3, 1]); hold on
end
plot(results_alt_1.continuous.x, results_prg_1.continuous.x, 'color', color_r, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_1.continuous.p, results_prg_1.continuous.p, 'color', color_r, 'linewidth', linewidth_cop, 'LineStyle', ':', 'Marker', 's', 'MarkerSize', markersize_cop, 'markerfacecolor', color_r, 'markeredgecolor', 'none', 'DisplayName', 'CoP');
plot(results_alt_2.continuous.x, results_prg_2.continuous.x, 'color', color_1, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_2.continuous.p, results_prg_2.continuous.p, 'color', color_1, 'linewidth', linewidth_cop, 'LineStyle', ':', 'Marker', 's', 'MarkerSize', markersize_cop, 'markerfacecolor', color_1, 'markeredgecolor', 'none', 'DisplayName', 'CoP');
plot(results_alt_3.continuous.x, results_prg_3.continuous.x, 'color', color_2, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
plot(results_alt_3.continuous.p, results_prg_3.continuous.p, 'color', color_2, 'linewidth', linewidth_cop, 'LineStyle', ':', 'Marker', 's', 'MarkerSize', markersize_cop, 'markerfacecolor', color_2, 'markeredgecolor', 'none', 'DisplayName', 'CoP');
xlabel('lateral position', 'fontsize', 14)
ylabel('anterior-posterior position', 'fontsize', 14)
% axis equal
xlim([-0.3 0.7])
ylim([-0.5 4])
if individual_figures && save_figures
    filename = ['fig_example_path' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(gcf)
end


if ~individual_figures && save_figures
    % save with labels
    filename_with = 'fig_trajectoryExample';
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename_with], '-djpeg', '-r300')
end

