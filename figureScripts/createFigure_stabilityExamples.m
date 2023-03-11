% flags
save_figure     = 0;

% set parameters
dt = 1e-4;
metronome_frequency = 1;
T_total_ref = 2;
T_total_pert = 8;
T_step = 1/(2*metronome_frequency);
z_c = 0.814;            % repeats the values in the simulation
g = 9.81;
omega = sqrt(g/z_c);
c = cosh(omega * T_step * 0.5);
s = sinh(omega * T_step * 0.5);

% define control parameters
b_o = 0.02;
b_p = 2.5;
b_d = 0.6;
% delta_b_d = [-0.2 -0.1 0 0.1 0.2];
delta_b_d = [-0.2 0 0.3 0.6];

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
gray = [1 1 1]*0.6;
colors = copper(number_of_systems);
markersize_cop = 12;
linewidth_cop = 2;
linewidth_com = 6;
fig = figure('position', [0, 0, 600 800]);
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'padding', 'compact');

% medial-lateral, position
nexttile(1); hold on
plot(results_ref.continuous.time, results_ref.continuous.x, 'color', gray, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM', 'displayname', 'reference');
% plot(results_ref.continuous.time, results_ref.continuous.p, 'color', gray, 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
for i_system = 1 : number_of_systems
    b_d_here = b_d + delta_b_d(i_system);
    plot(results{i_system}.continuous.time + T_total_ref, results{i_system}.continuous.x, 'color', colors(i_system, :), 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM', 'displayname', num2str(b_d_here));
%     plot(results{i_system}.continuous.time, results{i_system}.continuous.p, 'color', colors(i_system, :), 'linewidth', linewidth_cop, 'LineStyle', ':', 'DisplayName', 'CoP');
end
legend('location', 'best')
ylim([-0.2 0.4])
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral position', 'fontsize', 14)

% medial-lateral, velocity
nexttile(2); hold on
plot(results_ref.continuous.time, results_ref.continuous.v, 'color', gray, 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
for i_system = 1 : number_of_systems
    plot(results{i_system}.continuous.time + T_total_ref, results{i_system}.continuous.v, 'color', colors(i_system, :), 'linewidth', linewidth_com, 'LineStyle', '-', 'DisplayName', 'CoM');
end
ylim([-2 2])
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral velocity', 'fontsize', 14)



if save_figure
    % save with labels
    filename_with = 'fig_stabilityExamples';
    print(gcf, ['..' filesep 'figures' filesep filename_with], '-djpeg', '-r300')
%     filename_without = 'figX_example';

%     % remove text and marks to save data lines only
%     set(get(ax, 'xaxis'), 'visible', 'off');
%     set(get(ax, 'yaxis'), 'visible', 'off');
%     set(get(ax, 'xlabel'), 'visible', 'off');
%     set(get(ax, 'ylabel'), 'visible', 'off');
%     set(get(ax, 'title'), 'visible', 'off');
%     set(ax, 'xticklabel', '');
%     set(ax, 'yticklabel', '');
%     set(ax, 'position', [0 0 1 1]);
%     legend(ax, 'hide');
%     print(gcf, filename_without, '-djpeg', '-r300')
end

