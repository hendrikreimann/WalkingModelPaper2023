function ...
drawSlipPhasePortrait ...
  ( ...
    ax, ...
    p, ...
    x_min, ...
    x_max, ...
    v_min, ...
    v_max, ...
    orbits_intercepts_x, ...
    orbits_intercepts_v, ...
    omega, ...
    graphic_specs, ...
    arrow_specs ...
  )

    % default inputs
    if nargin < 10
        graphic_specs = struct;
        graphic_specs.linewidth = 0.5;
        graphic_specs.orbit_resolution = 100;
        graphic_specs.line_color = lines(1);
        graphic_specs.diagonal_color = lines(1);
        graphic_specs.diagonal_width = 2;
    end
    if nargin < 11
        arrow_specs = struct;
        arrow_specs.dt = 0.05;
        arrow_specs.width = 8;
        arrow_specs.length = 10;
        arrow_specs.color = lines(1);
    end

    % draw diagonals
    if -omega*(x_min-p) > v_max
        upper_left_x = -(1/omega * v_max - p);
        upper_left_v = v_max;
    else
        upper_left_x = x_min;
        upper_left_v = -(upper_left_x-p) * omega;
    end
    
    if -omega*(x_max-p) < v_min
        lower_right_x = -1/omega * v_min + p;
        lower_right_v = v_min;
    else
        lower_right_x = x_max;
        lower_right_v = -(lower_right_x-p) * omega;
    end
    plot([upper_left_x lower_right_x], [upper_left_v lower_right_v], 'linewidth', graphic_specs.diagonal_width, 'color', graphic_specs.diagonal_color)
    
    if omega*(x_max-p) > v_max
        upper_right_x = 1/omega * v_max + p;
        upper_right_v = v_max;
    else
        upper_right_x = x_max;
        upper_right_v = (upper_right_x-p) * omega;
    end
    
    if omega*(x_min-p) < v_min
        lower_left_x = 1/omega * v_min + p;
        lower_left_v = v_min;
    else
        lower_left_x = x_min;
        lower_left_v = (lower_left_x-p) * omega;
    end
    plot([lower_left_x upper_right_x], [lower_left_v upper_right_v], 'linewidth', graphic_specs.diagonal_width, 'color', graphic_specs.diagonal_color)
    
    % create orbits intersecting x-axis
    for i_orbit = 1 : length(orbits_intercepts_x)
        x_0 = orbits_intercepts_x(i_orbit);
        v_0 = 0;
        
        % find starting time
        if (x_0-p) < 0
            xlim_here = x_min;
            t_min_x = -acosh((xlim_here-p)/(x_0-p)) / omega;
            t_max_x = acosh((xlim_here-p)/(x_0-p)) / omega;
            t_min_v = asinh(v_max/(omega*(x_0-p))) / omega;
            t_max_v = asinh(v_min/(omega*(x_0-p))) / omega;

            t_min = max(t_min_x, t_min_v);
            t_max = min(t_max_x, t_max_v);
        elseif (x_0-p) > 0
            xlim_here = x_max;
            t_min_x = -acosh((xlim_here-p)/(x_0-p)) / omega;
            t_max_x = acosh((xlim_here-p)/(x_0-p)) / omega;
            t_min_v = asinh(v_min/(omega*(x_0-p))) / omega;
            t_max_v = asinh(v_max/(omega*(x_0-p))) / omega;

            t_min = max(t_min_x, t_min_v);
            t_max = min(t_max_x, t_max_v);
        end
        
        % plot
        time_here = linspace(t_min, t_max, graphic_specs.orbit_resolution);
        x_here = p + (x_0 - p)*cosh(omega*time_here) + v_0/omega*sinh(omega*time_here);
        v_here = (x_0 - p)*sinh(omega*time_here)*omega + v_0*cosh(omega*time_here);
        plot(ax, x_here, v_here, 'linewidth', graphic_specs.linewidth, 'color', graphic_specs.line_color);
        
        % arrows
        arrow_times_pos = arrow_specs.dt*0.5 : arrow_specs.dt : t_max;
        arrow_times_neg = -arrow_specs.dt*0.5 : -arrow_specs.dt : t_min;
        arrow_times = [flip(arrow_times_neg) arrow_times_pos];
        anchors_x = p + (x_0 - p)*cosh(omega*arrow_times) + v_0/omega*sinh(omega*arrow_times);
        anchors_v = (x_0 - p)*sinh(omega*arrow_times)*omega + v_0*cosh(omega*arrow_times);
        plot(ax, anchors_x, anchors_v, 'o', 'markersize', arrow_specs.width, 'markerfacecolor', arrow_specs.color, 'color', 'none');
    end

    % create orbits intersecting v-axis
    for i_orbit = 1 : length(orbits_intercepts_v)
        x_0 = p;
        v_0 = orbits_intercepts_v(i_orbit);
        
        % find starting time
        if v_0 < 0
            vlim_here = v_min;
            t_min_x = -acosh(vlim_here/v_0) / omega;
            t_max_x = acosh(vlim_here/v_0) / omega;
            t_min_v = asinh((omega*(x_max-p))/v_0) / omega;
            t_max_v = asinh((omega*(x_min-p))/v_0) / omega;
            
            t_min = max(t_min_x, t_min_v);
            t_max = min(t_max_x, t_max_v);
        elseif v_0 > 0
            vlim_here = v_max;
            t_min_x = -acosh(vlim_here/v_0) / omega;
            t_max_x = acosh(vlim_here/v_0) / omega;
            t_min_v = asinh((omega*(x_min-p))/v_0) / omega;
            t_max_v = asinh((omega*(x_max-p))/v_0) / omega;
            
            t_min = max(t_min_x, t_min_v);
            t_max = min(t_max_x, t_max_v);
        end
                
        % plot
        time_here = linspace(t_min, t_max, graphic_specs.orbit_resolution);
        x_here = p + (x_0 - p)*cosh(omega*time_here) + v_0/omega*sinh(omega*time_here);
        v_here = (x_0 - p)*sinh(omega*time_here)*omega + v_0*cosh(omega*time_here);
        plot(ax, x_here, v_here, 'linewidth', graphic_specs.linewidth, 'color', graphic_specs.line_color);

        % arrows
        arrow_times_pos = arrow_specs.dt*0.5 : arrow_specs.dt : t_max;
        arrow_times_neg = -arrow_specs.dt*0.5 : -arrow_specs.dt : t_min;
        arrow_times = [flip(arrow_times_neg) arrow_times_pos];
        anchors_x = p + (x_0 - p)*cosh(omega*arrow_times) + v_0/omega*sinh(omega*arrow_times);
        anchors_v = (x_0 - p)*sinh(omega*arrow_times)*omega + v_0*cosh(omega*arrow_times);
        
        plot(ax, anchors_x, anchors_v, 'o', 'markersize', arrow_specs.width, 'markerfacecolor', arrow_specs.color, 'color', 'none');
    end

    


