% flags
individual_figures          = 1;
save_figures                = 1;
labels                      = 'on';

% suppress warnings from flags
%#ok<*UNRCH> 

% visualization parameters
figure_width = 600;
figure_height = 300;

linewidth = 8;

color_r = [1 1 1]*0.6;
color_1 = hex2rgb('#20B309'); % green
color_2 = hex2rgb('#00A1E6'); % blue
color_3 = hex2rgb('#FFC51F'); % yellow
color_4 = hex2rgb('#CC3D0A'); % red
color_5 = hex2rgb('#880ACC'); % purple
colors = [color_1; color_2; color_3; color_5];

% set parameters
dt = 1e-4;
metronome_frequency = 1;
T_total_ref = 1;
T_total_pert = 9;
T_step = 1/(2*metronome_frequency);
z_c = 1;
g = 9.81;
omega = sqrt(g/z_c);
c = cosh(omega * T_step * 0.5);
s = sinh(omega * T_step * 0.5);

% define control parameters
b_o = 0.02;
b_p = 2.5;
b_d = 0.6;
delta_b_d = [-0.3 0 0.5 0.9];

% define initial conditions
p_0 = 0.0;
x_ref = b_o / (2*c - b_p);
v_ref = 0;
delta_x_0 = 0.01;
delta_v_0 = 0.0;

% define control structs
parameters.control = struct;
parameters.control.T_step = T_step;
parameters.control.b_o = b_o;
parameters.control.b_p = b_p;
parameters.control.b_d = b_d;
parameters.biomechanics = struct;
parameters.biomechanics.z_c = z_c;
parameters.biomechanics.g = g;

% solve
results_ref = simulateSlipWalker(T_total_ref, parameters, x_ref, v_ref, p_0, dt, 'alternating');
number_of_systems = numel(delta_b_d);
results = cell(number_of_systems, 1);
for i_system = 1 : number_of_systems
    this_delta_b_d = delta_b_d(i_system);
    parameters.control.b_d = b_d + this_delta_b_d;
    results{i_system} = simulateSlipWalker(T_total_pert, parameters, x_ref + delta_x_0, v_ref + delta_v_0, p_0, dt, 'alternating');
end

% visualize
if ~individual_figures
    fig = figure('position', [0, 0, 600 800]);
    tiledlayout(2, 1, 'TileSpacing', 'Compact', 'padding', 'compact');
end

% position
if individual_figures
    figure('position', [100 figure_height+180 figure_width figure_height]); hold on
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

plot(results_ref.continuous.time, results_ref.continuous.x, 'color', color_r, 'linewidth', linewidth, 'LineStyle', '-', 'DisplayName', 'CoM', 'displayname', 'reference');
for i_system = 1 : number_of_systems
    b_d_here = b_d + delta_b_d(i_system);
    plot(results{i_system}.continuous.time + T_total_ref, results{i_system}.continuous.x, 'color', colors(i_system, :), 'linewidth', linewidth, 'LineStyle', '-', 'DisplayName', 'CoM', 'displayname', num2str(b_d_here));
end
ylim([-0.1 0.5])
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral position', 'fontsize', 14)
if strcmp(labels, 'on')
    legend('location', 'best')
end

if individual_figures && save_figures
    filename = ['fig_stabilityExample_pos_' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
    close(gcf)
end

% velocity
if individual_figures
    figure('position', [100 100 figure_width figure_height]); hold on
    if strcmp(labels, 'on')
        label = 'withLabels';
        legend('location', 'best')
    elseif strcmp(labels, 'off')
        set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        label = 'noLabels';
    end
else
    nexttile(2); hold on
end

plot(results_ref.continuous.time, results_ref.continuous.v, 'color', color_r, 'linewidth', linewidth, 'LineStyle', '-', 'DisplayName', 'CoM', 'displayname', 'reference');
for i_system = 1 : number_of_systems
    b_d_here = b_d + delta_b_d(i_system);
    plot(results{i_system}.continuous.time + T_total_ref, results{i_system}.continuous.v, 'color', colors(i_system, :), 'linewidth', linewidth, 'LineStyle', '-', 'DisplayName', 'CoM vel', 'displayname', num2str(b_d_here));
end
ylim([-2 2])
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral position', 'fontsize', 14)

if individual_figures && save_figures
    filename = ['fig_stabilityExample_vel_' label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
    close(gcf)
end



