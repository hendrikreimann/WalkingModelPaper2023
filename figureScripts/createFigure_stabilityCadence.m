% flags
save_figure             = 0;
normalize_by_leg_length = 1;
report_r_square         = 0;
labels                  = 'on';

% suppress warnings from flags
%#ok<*UNRCH> 

% parameters
cadences_to_use = [80 110];
b_p_min = 0;
b_p_max = 5;
b_v_min = 0;
b_v_max = 1.5;
z_c = 0.77;
g = 9.81;

% colors
color_80 = hex2rgb('#009ADE');
color_110 = hex2rgb('#FFC61E');
% color_2 = hex2rgb('#00CD6C');
% color_3 = hex2rgb('#AF58BA');

colors = [color_80; color_110;];

% dependent variables
minute_to_second = 1/60;
number_of_cadences = length(cadences_to_use);


% load and extract data
data_file = '../../data/foot_placement_ml_by_com_50_CAD.csv';
human_data = readtable(data_file);
beta_p_data = cell(number_of_cadences, 1);
beta_v_data = cell(number_of_cadences, 1);
r_square_left_data = cell(number_of_cadences, 1);
r_square_right_data = cell(number_of_cadences, 1);
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
    r_square_left = table2array(human_data(indicator_left, 'R_square'));
    r_square_right = table2array(human_data(indicator_right, 'R_square'));
    
    if isequal(subjects_left, subjects_right)
        beta_p_left = beta_p_left_unsorted;
        beta_p_right = beta_p_right_unsorted;
        beta_v_left = beta_v_left_unsorted;
        beta_v_right = beta_v_right_unsorted;
        subjects = subjects_left;
    else
        error('Subject order for left and right not the same.')
    end
    beta_p = (beta_p_left + beta_p_right) * 0.5;
    beta_v = (beta_v_left + beta_v_right) * 0.5;

    beta_p_data{i_cadence} = beta_p;
    beta_v_data{i_cadence} = beta_v;

    r_square_left_data{i_cadence} = r_square_left;
    r_square_right_data{i_cadence} = r_square_right;
    
end

% report r square
if report_r_square
    r_square_left_80_min = min(r_square_left_data{1})
    r_square_left_110_min = min(r_square_left_data{2})
    r_square_left_80_max = max(r_square_left_data{1})
    r_square_left_110_max = max(r_square_left_data{2})

    r_square_right_80_min = min(r_square_right_data{1})
    r_square_right_110_min = min(r_square_right_data{2})
    r_square_right_80_max = max(r_square_right_data{1})
    r_square_right_110_max = max(r_square_right_data{2})

    r_square_80_min = min([r_square_left_data{1}; r_square_right_data{1}])
    r_square_110_min = min([r_square_left_data{2}; r_square_right_data{2}])
    r_square_80_max = max([r_square_left_data{1}; r_square_right_data{1}])
    r_square_110_max = max([r_square_left_data{2}; r_square_right_data{2}])

    return
end

% normalize
if normalize_by_leg_length
    data_file = '../../data/legLength_CAD.mat';
    loaded_data = load(data_file);
    length_data = loaded_data.data_table;
    
    for i_subject = 1 : length(subjects)
        this_subject = subjects{i_subject};
        this_subject_leg_length = length_data{strcmp(length_data.subject, this_subject), "leg_length"};
        for i_cadence = 1 : number_of_cadences
            beta_p_data{i_cadence}(i_subject) = beta_p_data{i_cadence}(i_subject) * 1/this_subject_leg_length;
            beta_v_data{i_cadence}(i_subject) = beta_v_data{i_cadence}(i_subject) * 1/sqrt(this_subject_leg_length);
        end
    end
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
patch_handles = [];
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
    h = patch([P1(1) P2(1) P3(1)], [P1(2) P2(2) P3(2)], [0 0 0], 'facecolor', lightenColor(this_color, 0.5), 'edgecolor', 'none');
    patch_handles = [patch_handles, h]; %#ok<AGROW> 

    % plot data
    markersize = 10;
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

for i_cadence = 1 : number_of_cadences
    uistack(patch_handles(i_cadence), 'bottom')
end
% save
if save_figure
    % save with labels
    filename = ['fig_stabilityCadence_' file_label];
    print(gcf, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r600')
end
 






























