% flags
save_figures                = 0;
run_stats                   = 1;
labels                      = 'on';

% suppress warnings from flags
%#ok<*UNRCH> 

% load and extract data
data_file = '../../data/16-Nov-2022_R2_stepwidth.mat';
load(data_file)
soi=25;
beta_p = squeeze(beta1(:,soi,:));
beta_p_normal = beta_p(:,3);
beta_p_lesschuh = beta_p(:,4);

beta_d = squeeze(beta2(:,soi,:));
beta_d_normal = beta_d(:,3);
beta_d_lesschuh = beta_d(:,4);

% stats
if run_stats
    [hypothesis_p, p_value_p] = ttest(beta_p_normal, beta_p_lesschuh);
    [hypothesis_v, p_value_v] = ttest(beta_d_normal, beta_d_lesschuh);
end

% define visualization stuff
figure_width = 400;
figure_height = 400;
xlimits = [0.2 3.8];

color_r = [1 1 1]*0.6;
color_1 = hex2rgb('#20B309'); % green
color_2 = hex2rgb('#00A1E6'); % blue

marker = 'o';
markersize = 12;
jitter_width = 0.1;
jitter = (rand(size(beta_p_normal))-0.5) * jitter_width;
width = 1.0;

% pos
fig_pos = figure('position', [100 100 figure_width figure_height], 'renderer', 'painters'); hold on
if strcmp(labels, 'on')
    label = 'withLabels';
elseif strcmp(labels, 'off')
    set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
    set(get(gca, 'xaxis'), 'visible', 'off');
    set(get(gca, 'yaxis'), 'visible', 'off');
    label = 'noLabels';
end
xlim(xlimits)
bar(1:3,[mean(beta_p_normal) NaN mean(beta_p_lesschuh)], width, 'linewidth', 2)
for i_line = 1 : length(beta_p_normal)
    plot([1 3] + jitter(i_line), [beta_p_normal(i_line, :) beta_p_lesschuh(i_line, :)],'Color', color_r, 'linewidth', 1)
end
plot(1+jitter, beta_p_normal, 'marker', marker, 'linestyle', 'none', 'MarkerFaceColor', color_r, 'color', 'k', 'MarkerSize', markersize, 'linewidth', 2)
plot(3+jitter, beta_p_lesschuh, 'marker', marker, 'linestyle', 'none', 'MarkerFaceColor', color_r, 'color', 'k', 'MarkerSize', markersize, 'linewidth', 2)

xticks(1:3)
xticklabels({'Normal',[],'Lesschuh'})
h=gca; h.XAxis.TickLength = [0 0];
box off

title(sprintf('\\beta_p at %.0f%% of step',soi/50*100))
ylabel('\beta_1')

if save_figures
    filename = ['fig_beta_p_' label];
    print(fig_pos, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(fig_pos)
end

% vel
fig_vel = figure('position', [100 100 figure_width figure_height], 'renderer', 'painters'); hold on
if strcmp(labels, 'on')
    label = 'withLabels';
elseif strcmp(labels, 'off')
    set(gca, 'position', [0 0 1 1], 'xtick', [], 'ytick', []);
    set(get(gca, 'xaxis'), 'visible', 'off');
    set(get(gca, 'yaxis'), 'visible', 'off');
    label = 'noLabels';
end
xlim(xlimits)
bar(1:3,[mean(beta_d_normal) NaN mean(beta_d_lesschuh)], width, 'linewidth', 2)
for i_line = 1 : length(beta_d_normal)
    plot([1 3] + jitter(i_line), [beta_d_normal(i_line, :) beta_d_lesschuh(i_line, :)],'Color', color_r, 'linewidth', 1)
end
plot(1+jitter, beta_d_normal, 'marker', marker, 'linestyle', 'none', 'MarkerFaceColor', color_r, 'color', 'k', 'MarkerSize', markersize, 'linewidth', 2)
plot(3+jitter, beta_d_lesschuh, 'marker', marker, 'linestyle', 'none', 'MarkerFaceColor', color_r, 'color', 'k', 'MarkerSize', markersize, 'linewidth', 2)

xticks(1:3)
xticklabels({'Normal',[],'Lesschuh'})
h=gca; h.XAxis.TickLength = [0 0];
box off

title(sprintf('\\beta_d at %.0f%% of step',soi/50*100))
ylabel('\beta_1')

if save_figures
    filename = ['fig_beta_d_' label];
    print(fig_vel, ['..' filesep '..' filesep 'figures_raw' filesep filename], '-djpeg', '-r300')
    close(fig_vel)
end




