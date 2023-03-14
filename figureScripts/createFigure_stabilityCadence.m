
% flags
save_figure             = 1;
labels                  = 'on';

% parameters
cadences_to_use = [80 110];
b_p_min = 0;
b_p_max = 5;
b_v_min = 0;
b_v_max = 1.5;
z_c = 0.814;
g = 9.81;

% colors
color_80 = hex2rgb('#009ADE');
color_110 = hex2rgb('#FFC61E');
% color_2 = hex2rgb('#00CD6C');
% color_3 = hex2rgb('#AF58BA');

colors = [color_80; color_110; color_2; color_3];

% dependent variables
minute_to_second = 1/60;
number_of_cadences = length(cadences_to_use);


% load and extract data
data_file = '../../data/foot_placement_ml_by_com_50_CAD.csv';
human_data = readtable(data_file);
beta_p_data = cell(number_of_cadences, 1);
beta_v_data = cell(number_of_cadences, 1);
for i_cadence = 1 : number_of_cadences
    this_cadence = cadences_to_use(i_cadence);
    
    indicator_left = strcmp(table2cell(human_data(:, 'cadence')), [num2str(this_cadence) 'BPM']) & strcmp(table2cell(human_data(:, 'direction')), 'STIM_TOWARDS') & strcmp(table2cell(human_data(:, 'trigger_foot')), 'TRIGGER_LEFT');
    indicator_right = strcmp(table2cell(human_data(:, 'cadence')), [num2str(this_cadence) 'BPM']) & strcmp(table2cell(human_data(:, 'direction')), 'STIM_TOWARDS') & strcmp(table2cell(human_data(:, 'trigger_foot')), 'TRIGGER_RIGHT');
    beta_p_left_unsorted = table2array(human_data(indicator_left, 'foot_placement_ml_by_com_50_slope_1'));
    beta_p_right_unsorted = table2array(human_data(indicator_right, 'foot_placement_ml_by_com_50_slope_1'));
    beta_v_left_unsorted = table2array(human_data(indicator_left, 'foot_placement_ml_by_com_50_slope_2'));
    beta_v_right_unsorted = table2array(human_data(indicator_right, 'foot_placement_ml_by_com_50_slope_2'));
    subjects_left = table2cell(human_data(indicator_left, 'subject'));
    subjects_right = table2cell(human_data(indicator_right, 'subject'));
    
    if isequal(subjects_left, subjects_right)
        beta_p_left = beta_p_left_unsorted;
        beta_p_right = beta_p_right_unsorted;
        beta_v_left = beta_v_left_unsorted;
        beta_v_right = beta_v_right_unsorted;
    else
        error('Subject order for left and right not the same. Need to sort by subject, but this is not implemented yet.')
    end
    beta_p = (beta_p_left + beta_p_right) * 0.5;
    beta_v = (beta_v_left + beta_v_right) * 0.5;

    beta_p_data{i_cadence} = beta_p;
    beta_v_data{i_cadence} = beta_v;
end

% create figure
figure; ax = axes; hold on
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

% stability regions
for i_cadence = 1 : number_of_cadences
    % calculate dependent variables
    bpm = cadences_to_use(i_cadence);
    metronome_frequency = bpm * 0.5 * minute_to_second; % two beats per metronome cycle
    stride_time = metronome_frequency.^(-1);
    step_time = stride_time * 0.5;
    T_step = 1/(2*metronome_frequency);
    omega = sqrt(g/z_c);
    c = cosh(omega * T_step * 0.5);
    s = sinh(omega * T_step * 0.5);

    % color
    this_color = colors(i_cadence, :);

    % plot triangle
    P1 = [2*c, 2*s/omega];
    P2 = [2*s^2/c, 2*s/omega];
    P3 = [2*c, 2*c^2/(omega*s)];
    patch([P1(1) P2(1) P3(1)], [P1(2) P2(2) P3(2)], [0 0 0], 'facecolor', lightenColor(this_color, 0.5), 'edgecolor', 'none')

    % plot data
    markersize = 8;
    plot ...
      ( ...
        beta_p_data{i_cadence}, beta_v_data{i_cadence}, ...
        'o', ...
        'linewidth', 0.5, ...
        'MarkerEdgeColor', [1 1 1]*0.5, ...
        'MarkerFaceColor', colors(i_cadence, :), ...
        'markersize', markersize ...
      );

end

% save
if save_figure
    % save with labels
    filename = ['fig_stabilityCadence_' file_label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
end
 






























