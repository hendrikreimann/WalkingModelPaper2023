
angles = [0 2*pi/12 2*pi/6];
% angles = 2*pi/6;

targets = [0.3; 0.1];
% angles = atan2(target(2), target(1));

color = lines(1);
arrow_length = 50;
arrow_width = 5;


figure('Position', [0 0 1200 800]);
tiledlayout(3, 3);
ax1 = nexttile(1); hold on; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
ax2 = nexttile(2, [1 2]); hold on; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
ax3 = nexttile(4, [2 1]); hold on; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
ax4 = nexttile(5, [2 2]); hold on; xlim([-1.1 1.1]); ylim([-1.1 1.1]);

plot(ax1, cos(0 : 0.001 : 2*pi), sin(0 : 0.001 : 2*pi));
plot(ax2, cos(0 : 0.001 : 2*pi), sin(0 : 0.001 : 2*pi));
plot(ax3, cos(0 : 0.001 : 2*pi), sin(0 : 0.001 : 2*pi));
plot(ax4, cos(0 : 0.001 : 2*pi), sin(0 : 0.001 : 2*pi));

for i_angle = 1 : size(targets, 2)
    
    plotArrowhead(ax1, [0; 0], targets(:, i_angle), arrow_length, arrow_width, color)
    plotArrowhead(ax2, [0; 0], targets(:, i_angle), arrow_length, arrow_width, color)
    plotArrowhead(ax3, [0; 0], targets(:, i_angle), arrow_length, arrow_width, color)
    plotArrowhead(ax4, [0; 0], targets(:, i_angle), arrow_length, arrow_width, color)
    
    plot(ax1, [0 targets(1, i_angle)], [0, targets(2, i_angle)])
    plot(ax2, [0 targets(1, i_angle)], [0, targets(2, i_angle)])
    plot(ax3, [0 targets(1, i_angle)], [0, targets(2, i_angle)])
    plot(ax4, [0 targets(1, i_angle)], [0, targets(2, i_angle)])
end
