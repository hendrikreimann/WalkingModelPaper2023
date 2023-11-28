% flags
save_figures            = 0;
labels                  = 'on';
run_stats               = 1;
normalize_by_leg_length = 1;

% suppress warnings from flags
%#ok<*UNRCH> 

% parameters
cadences = 70 : 120;
z_c = 0.77;
g = 9.81;

% play around with different values of omega

% visualization parameters
figure_width = 600;
figure_height = 300;
gray = [0.4 0.4 0.4];
color_80 = hex2rgb('#009ADE');
color_110 = hex2rgb('#FFC61E');
% color_2 = hex2rgb('#00CD6C');
% color_3 = hex2rgb('#AF58BA');

colors = [color_80; color_110];

% dependent variables
minute_to_second = 1/60;
number_of_cadences = length(cadences);

% model
b_p_min = zeros(1, number_of_cadences);
b_p_max = zeros(1, number_of_cadences);
b_d_min = zeros(1, number_of_cadences);
b_d_max = zeros(1, number_of_cadences);
for i_cadence = 1 : number_of_cadences
    % calculate dependent variables
    bpm = cadences(i_cadence);
    metronome_frequency = bpm * 0.5 * minute_to_second; % two beats per metronome cycle
    stride_time = metronome_frequency.^(-1);
    T_step = 1/(2*metronome_frequency);
    omega = sqrt(g/z_c);
    c = cosh(omega * T_step * 0.5);
    s = sinh(omega * T_step * 0.5);

    % stability region vertices
    P1 = [2*c, 2*s/omega];
    P2 = [2*s^2/c, 2*s/omega];
    P3 = [2*c, 2*c^2/(omega*s)];

    b_p_min(i_cadence) = 2*s^2/c;
    b_p_max(i_cadence) = 2*c;
    b_d_min(i_cadence) = 2*s/omega;
    b_d_max(i_cadence) = 2*c^2/(omega*s);
end

% load and extract data
cadences_human = [80, 110];
number_of_cadences_human = length(cadences_human);
data_file = '../../data/foot_placement_ml_by_com_50_CAD.csv';
human_data = readtable(data_file);
beta_p_data = cell(number_of_cadences_human, 1);
beta_v_data = cell(number_of_cadences_human, 1);
cadence_actual_data = cell(number_of_cadences_human, 1);
for i_cadence = 1 : number_of_cadences_human
    this_cadence = cadences_human(i_cadence);
    
    indicator_left = strcmp(table2cell(human_data(:, 'cadence')), [num2str(this_cadence) 'BPM']) & strcmp(table2cell(human_data(:, 'direction')), 'STIM_TOWARDS') & strcmp(table2cell(human_data(:, 'trigger_foot')), 'TRIGGER_LEFT');
    indicator_right = strcmp(table2cell(human_data(:, 'cadence')), [num2str(this_cadence) 'BPM']) & strcmp(table2cell(human_data(:, 'direction')), 'STIM_TOWARDS') & strcmp(table2cell(human_data(:, 'trigger_foot')), 'TRIGGER_RIGHT');
    beta_p_left_unsorted = table2array(human_data(indicator_left, 'foot_placement_ml_by_com_50_slope_1'));
    beta_p_right_unsorted = table2array(human_data(indicator_right, 'foot_placement_ml_by_com_50_slope_1'));
    beta_v_left_unsorted = table2array(human_data(indicator_left, 'foot_placement_ml_by_com_50_slope_2'));
    beta_v_right_unsorted = table2array(human_data(indicator_right, 'foot_placement_ml_by_com_50_slope_2'));
    cadence_actual_left_unsorted = table2array(human_data(indicator_left, 'cadence_actual'));
    cadence_actual_right_unsorted = table2array(human_data(indicator_right, 'cadence_actual'));
    subjects_left = table2cell(human_data(indicator_left, 'subject'));
    subjects_right = table2cell(human_data(indicator_right, 'subject'));
    
    if isequal(subjects_left, subjects_right)
        beta_p_left = beta_p_left_unsorted;
        beta_p_right = beta_p_right_unsorted;
        beta_v_left = beta_v_left_unsorted;
        beta_v_right = beta_v_right_unsorted;
        cadence_actual_left = cadence_actual_left_unsorted;
        cadence_actual_right = cadence_actual_right_unsorted;
    else
        error('Subject order for left and right not the same. Need to sort by subject, but this is not implemented yet.')
    end
    beta_p = (beta_p_left + beta_p_right) * 0.5;
    beta_v = (beta_v_left + beta_v_right) * 0.5;
    cadence_actual = (cadence_actual_left + cadence_actual_right) * 0.5;

    beta_p_data{i_cadence} = beta_p;
    beta_v_data{i_cadence} = beta_v;
    cadence_actual_data{i_cadence} = cadence_actual;
end

% normalize
if normalize_by_leg_length
    data_file = '../../data/legLength_CAD.mat';
    loaded_data = load(data_file);
    length_data = loaded_data.data_table;
    
    for i_subject = 1 : length(subjects)
        this_subject = subjects{i_subject};
        this_subject_leg_length = length_data{strcmp(length_data.subject, this_subject), "leg_length"};
        for i_cadence = 1 : number_of_cadences_human
            beta_p_data{i_cadence}(i_subject) = beta_p_data{i_cadence}(i_subject) * 1/this_subject_leg_length;
            beta_v_data{i_cadence}(i_subject) = beta_v_data{i_cadence}(i_subject) * 1/sqrt(this_subject_leg_length);
        end
    end
end

% stats
if run_stats
    [hypothesis_p, p_value_p] = ttest(beta_p_data{1}, beta_p_data{2});
    [hypothesis_v, p_value_v] = ttest(beta_v_data{1}, beta_v_data{2});
end

% figure for b_p
fig_pos = figure('position', [100 figure_height+180 figure_width figure_height]); axes_pos = axes; hold on
if strcmp(labels, 'on')
    label = 'withLabels';
elseif strcmp(labels, 'off')
    set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
    set(get(gca, 'xaxis'), 'visible', 'off');
    set(get(gca, 'yaxis'), 'visible', 'off');
    label = 'noLabels';
end
ylim([0, 5]);
plot_handles = shadedErrorBar ...
  ( ...
    cadences, ...
    (b_p_min+b_p_max)*0.5, ...
    (b_p_max-b_p_min)*0.5, ...
    { ...
      'color', gray, ...
      'linewidth', 6 ...
    }, ...
    1, ...
    axes_pos ...
  );
delete(plot_handles.mainLine);
set(plot_handles.edge, 'color', 'none');

% figure for b_d
fig_vel = figure('position', [100 100 figure_width figure_height]); axes_vel = axes; hold on
if strcmp(labels, 'on')
    label = 'withLabels';
elseif strcmp(labels, 'off')
    set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
    set(get(gca, 'xaxis'), 'visible', 'off');
    set(get(gca, 'yaxis'), 'visible', 'off');
    label = 'noLabels';
end
ylim([0, 1.5]);
plot_handles = shadedErrorBar ...
  ( ...
    cadences, ...
    (b_d_min+b_d_max)*0.5, ...
    (b_d_max-b_d_min)*0.5, ...
    { ...
      'color', gray, ...
      'linewidth', 6 ...
    }, ...
    1, ...
    axes_vel ...
  );
delete(plot_handles.mainLine);
set(plot_handles.edge, 'color', 'none');


for i_cadence = 1 : number_of_cadences_human
    this_cadence_actual_data = cadence_actual_data{i_cadence};
    this_beta_p_data = beta_p_data{i_cadence};
    this_beta_v_data = beta_v_data{i_cadence};
    this_color = colors(i_cadence, :);

    markersize = 12;
    plot ...
      ( ...
        axes_pos, ...
        this_cadence_actual_data, this_beta_p_data, ...
        'o', ...
        'linewidth', 0.5, ...
        'MarkerEdgeColor', [1 1 1]*0.5, ...
        'MarkerFaceColor', colors(i_cadence, :), ...
        'markersize', markersize ...
      );
    
    plot ...
      ( ...
        axes_vel, ...
        this_cadence_actual_data, this_beta_v_data, ...
        'o', ...
        'linewidth', 0.5, ...
        'MarkerEdgeColor', [1 1 1]*0.5, ...
        'MarkerFaceColor', colors(i_cadence, :), ...
        'markersize', markersize ...
      );
    
end


if save_figures
    filename = ['fig_cadenceSweep_pos_' label];
    print(fig_pos, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
    close(fig_pos)

    filename = ['fig_cadenceSweep_vel_' label];
    print(fig_vel, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
    close(fig_vel)
end



