% simulate SLIP walker with linear step placement feedback

% flags
save_figure                 = 0;
calculate_parameter_surf    = 1;
plot_points_and_lines_on_surf = 1;


% ----------------------------------------------------------------------------------------------------------------------
% pasted from createFigure_stabilityExamples.m
% ----------------------------------------------------------------------------------------------------------------------
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
% ----------------------------------------------------------------------------------------------------------------------


% 
b_p_number_of_steps = 401;
b_v_number_of_steps = 401;
% 
metronome_frequency = 1;
% 
% metronome_frequency = 1.5;
b_p_min = 0;
b_p_max = 5;
b_v_min = 0;
b_v_max = 2;

% define noise
noise = struct;
noise.sigma_step = 0;
noise.sigma_sensor_p = 0.001;
noise.sigma_sensor_v = 0.001;
noise.sigma_time = 0;
noise.rngSeed = 0;
noise.stream = RandStream('mt19937ar');

% calculate numeric solutions for a given x_0
results = simulateSlipWalker(T_total_ref, parameters, x_ref, v_ref, p_0, dt);
time = results.continuous.time;
step_index = results.indices.step_indices_continuous(1);
time_before_step = time(1 : step_index-1);
time_after_step = time(step_index : end);
time_after_step_rel = time_after_step - time_after_step(1) + time(1);
x_m_numeric = results.continuous.x(end); % x_m_numeric should be equal to x_m. This assumes the T = T_step
v_m_numeric = results.continuous.v(end); % v_m_numeric should be equal to v_m. This assumes the T = T_step
p_m_numeric = results.continuous.p(end);
q_m_numeric = x_m_numeric - p_m_numeric;


M = [ ...
      (- c*b_p + c^2 + s^2),  - (c*b_d - 2/omega*c*s); ...
      (2*c - b_p)*omega*s, (- b_d*omega*s + c^2 + s^2) ...
    ];

if calculate_parameter_surf
    % sample the region of the parameter space where both Eigenvalues are <1 absolute
    
    b_p_grid = linspace(b_p_min, b_p_max, b_p_number_of_steps);
    b_v_grid = linspace(b_v_min, b_v_max, b_v_number_of_steps);
    steps_per_minute = metronome_frequency * 120;
    
    
    [B_P,B_V] = meshgrid(b_p_grid,b_v_grid);
    Z_rho = zeros(size(B_P));
    C_rho = zeros(size(B_P));
    alpha_data_rho = zeros(size(B_P));
    
    Z_D = zeros(size(B_P));
    C_D = zeros(size(B_P));
    
    for i_p = 1 : b_p_number_of_steps
        b_p_here = b_p_grid(i_p);
        for i_v = 1 : b_v_number_of_steps
            b_v_here = b_v_grid(i_v);
            M_here = ...
                [ ...
                  (- c*b_p_here + c^2 + s^2),  - (c*b_v_here - 2/omega*c*s); ...
                  (2*c - b_p_here)*omega*s, (- b_v_here*omega*s + c^2 + s^2) ...
                ];
            
            % Eigendecomposition to get spectral norm rho(M) (https://en.wikipedia.org/wiki/Convergent_matrix)
            D = eig(M_here);
            rho = max(abs(D));
            Z_rho(i_v, i_p) = rho;
            Z_D(i_v, i_p) = trace(M_here)^2 - 4*det(M_here);
            
            if rho < 1
                C_rho(i_v, i_p) = (1-rho) * 10;
            else
%                 C_rho(i_v, i_p) = 0;
                
                if Z_D(i_v, i_p) < 0
                    C_rho(i_v, i_p) = 0.3;
                else
                    C_rho(i_v, i_p) = 0;
                end
                
            end
            
            % alternative coloring scheme for poster
            if rho < 1
                C_rho(i_v, i_p) = (1-rho) * 10;
                C_rho(i_v, i_p) = 1-rho;
%                 alpha_data_rho(i_v, i_p) = ((1-rho) + 1) * 0.5;
            else
                C_rho(i_v, i_p) = 0;
                alpha_data_rho(i_v, i_p) = 1;
                
            end

            
            
            
            
            
            if Z_D(i_v, i_p) < 0
                C_D(i_v, i_p) = 1;
            else
                C_D(i_v, i_p) = 0;
            end
            
        end
    end
    
%     C = X.*Y;
    figure; ax = axes; hold on; title(['step time = ' num2str(parameters.control.T_step) ', Cadence = ' num2str(steps_per_minute) ' steps/min'])
    xlabel('$b_p$'); ylabel('$b_v$')
    surf(B_P, B_V, -Z_rho, C_rho, 'AlphaData', alpha_data_rho, 'facecolor', 'interp')
%     colormap(ax, turbo)
    shading interp
    colorbar
    
    % add delimiter lines
    gray = [0.5 0.5 0.5];
    linewidth = 3;
    plot(2*c*[1 1], [b_v_min, b_v_max], 'color', gray, 'linewidth', linewidth);
    plot([b_p_min, b_p_max], 2*s/omega*[1 1], 'color', gray, 'linewidth', linewidth);

    plot([b_p_min, b_p_max], [b_p_min, b_p_max]*c/(s*omega), 'color', gray, 'linewidth', linewidth);
    xlim([b_p_min, b_p_max]);
    ylim([b_v_min, b_v_max]);
    
%     figure; ax = axes; hold on; title(['step time = ' num2str(control.T_step) ', Cadence = ' num2str(steps_per_minute) ' steps/min'])
%     xlabel('$b_p$'); ylabel('$b_v$')
%     surf(B_P, B_V, Z_D, C_D)
%     colormap(ax, turbo)
%     shading interp
    if plot_points_and_lines_on_surf
        plot3(2*c, 2*s/omega, -1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [1 0.5 0])
        plot3(2*s^2/c, 2*s/omega, -1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 1 0])
        plot3(2*c, 2*c^2/(omega*s), -1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 1 0])

        plot3((s^2+c^2)/c, (s^2+c^2)/(omega*s), -1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [1 0.5 0])

        plot3([2*c (s^2+c^2)/c], [2*s/omega (s^2+c^2)/(omega*s)], [0 0], '-', 'MarkerEdgeColor', 'none', 'markerfacecolor', [1 0.5 0])

        mu = 0.5
        plot3(2*c + mu*((s^2+c^2)/c - 2*c), 2*s/omega + mu*((s^2+c^2)/(omega*s) - 2*s/omega), 0, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0.5 1])
        plot3(0.5*(4*c^2-1)/c, 0.5*(4*s^2+1)/(s*omega), 0, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 1 0])

        % test out quadratic parameterization
        u = 0;
        f_c = -u^2 + 2*u + 3;
        f_s = u^2 - 2*u + 1;
        g_c = u^2 + 2*u + 1;
        g_s = -u^2 - 2*u + 3;    
        P_p = (f_c*c^2 + f_s*s^2) / (2*c);
        P_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(P_p, P_d, 1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0 1])

        u = 1;
        f_c = -u^2 + 2*u + 3;
        f_s = u^2 - 2*u + 1;
        g_c = u^2 + 2*u + 1;
        g_s = -u^2 - 2*u + 3;
        P_p = (f_c*c^2 + f_s*s^2) / (2*c);
        P_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(P_p, P_d, 1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0 1])

        u = -1;
        f_c = -u^2 + 2*u + 3;
        f_s = u^2 - 2*u + 1;
        g_c = u^2 + 2*u + 1;
        g_s = -u^2 - 2*u + 3;
        P_p = (f_c*c^2 + f_s*s^2) / (2*c);
        P_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(P_p, P_d, 1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0 1])

        u = 0.5;
        f_c = -u^2 + 2*u + 3;
        f_s = u^2 - 2*u + 1;
        g_c = u^2 + 2*u + 1;
        g_s = -u^2 - 2*u + 3;
        P_p = (f_c*c^2 + f_s*s^2) / (2*c);
        P_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(P_p, P_d, 1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0 1])

        u = 0.2;
        f_c = -u^2 + 2*u + 3;
        f_s = u^2 - 2*u + 1;
        g_c = u^2 + 2*u + 1;
        g_s = -u^2 - 2*u + 3;
        P_p = (f_c*c^2 + f_s*s^2) / (2*c);
        P_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(P_p, P_d, 1, 'o', 'MarkerEdgeColor', 'none', 'markerfacecolor', [0 0 1])

        % parameterize line
        u = [-3 : 0.01 : 3];
        f_c = -u.^2 + 2*u + 3;
        f_s = u.^2 - 2*u + 1;
        g_c = u.^2 + 2*u + 1;
        g_s = -u.^2 - 2*u + 3;
        L_p = (f_c*c^2 + f_s*s^2) / (2*c);
        L_d = (g_c*c^2 + g_s*s^2) / (2*s*omega);
        plot3(L_p, L_d, ones(size(u))*2, '-', 'linewidth', 2)
    end
    
    if save_figure
        % save with labels
        filename_with = 'fig7_withLabels';
        print(gcf, filename_with, '-djpeg', '-r300')
        filename_without = 'fig7_noLabels';

        % remove text and marks to save data lines only
        set(get(gca, 'xaxis'), 'visible', 'off');
        set(get(gca, 'yaxis'), 'visible', 'off');
        set(get(gca, 'xlabel'), 'visible', 'off');
        set(get(gca, 'ylabel'), 'visible', 'off');
        set(get(gca, 'title'), 'visible', 'off');
        set(gca, 'xticklabel', '');
        set(gca, 'yticklabel', '');
        set(gca, 'position', [0 0 1 1]);
        legend(gca, 'hide');
        print(gcf, filename_without, '-djpeg', '-r300')
        
    end
    
end







