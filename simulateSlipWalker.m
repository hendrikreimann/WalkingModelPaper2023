function ...
      [ ...
        results, ...
        control ...
      ] = ...
    simulateSlipWalker ... 
      ( ...
        T_total, ...
        parameters, ...
        x_init, ...
        v_init, ...
        p_init, ...
        dt, ...
        mode ...
      )

    % simulate walker walking system with proportional derivative control for foot placement
    if nargin < 1
        T_total = 120;
    end
    
    % set parameters
    if nargin < 6
        dt = 1e-3;
    end
    
    if nargin < 7
        mode = 'progressive';
    end
    
    time = dt : dt : T_total;
    number_of_time_steps = length(time);
      
    parameters.control.metronome_frequency = 1 / (2*parameters.control.T_step);
    
    % initialize discrete variables (one entry per step)
    step_duration_n = [];
    step_time_n = [];
    foot_placement_n = [];
    step_distance_n = [];
    com_from_cop_at_step_end_n = [];
    com_vel_at_step_end_n = [];
    com_from_cop_at_midstance_n = [];
    com_vel_at_midstance_n = [];

    n = 0;
    
    % initialize time series for continuous variables
    x_timeseries = zeros(number_of_time_steps, 1);
    v_timeseries = zeros(number_of_time_steps, 1);
    a_timeseries = zeros(number_of_time_steps, 1);
    u_timeseries = zeros(number_of_time_steps, 1);
    q_timeseries = zeros(number_of_time_steps, 1);
    p_timeseries = zeros(number_of_time_steps, 1);
    xcom_timeseries = zeros(number_of_time_steps, 1);
    metronome_timeseries = cos(2 * pi * time * parameters.control.metronome_frequency);
    metronome_timeseries(metronome_timeseries>=0) = 1;
    metronome_timeseries(metronome_timeseries<0) = 0;
    metronome_reset_indices = find(diff(metronome_timeseries)) + 1;
    metronome_phase = wrapToPi(2 * pi * time * parameters.control.metronome_frequency);
    this_metronome_state = metronome_timeseries(1);
        
    % find midstance indices
    step_indices = find(diff(metronome_timeseries));
    if numel(step_indices) == 1
        midstance_offsets = step_indices;
    else
        midstance_offsets = round(diff(step_indices) * 0.5);
    end
    median_midstance_offset = median(midstance_offsets);
    midstance_indices = [1 step_indices + [midstance_offsets median_midstance_offset]];
    
    step_number_timeseries = zeros(number_of_time_steps, 1);
    step_number_timeseries(1) = n;

    % initialize counters for different sides
    step_indices_continuous = [];
    left_step_indices_continuous = [];
    right_step_indices_continuous = [];
    step_indices_discrete = [];
    left_step_indices_discrete = [];
    right_step_indices_discrete = [];

    % define system dynamics
    omega = sqrt(parameters.biomechanics.g/parameters.biomechanics.z_c);
    
    % initialize system
    p = p_init;
    x = x_init;
    v = v_init;
    
    % define control
    com_from_step_midstance = 0;
    com_vel_midstance = 0;
    for i_time = 1 : number_of_time_steps
        % time
        t = time(i_time);
        last_metronome_state = this_metronome_state;
        this_metronome_state = metronome_timeseries(i_time);

        % at midstance, update CoM kinematic state estimates
        if ismember(i_time, midstance_indices)
            % sense state variables
            com_from_step_midstance = x - p;
            com_vel_midstance = v;
            
            % store sensed variables
            com_from_cop_at_midstance_n = [com_from_cop_at_midstance_n, com_from_step_midstance]; %#ok<AGROW>
            com_vel_at_midstance_n = [com_vel_at_midstance_n, com_vel_midstance]; %#ok<AGROW>
        end
        
        % step control
        if last_metronome_state ~= this_metronome_state
            % use CoM state at midstance as sensory input for control
            if strcmp(mode, 'progressive')
                % Equation 3 (eq_footPlacement_ap)
                p_new = p ...
                        + parameters.control.b_o ...
                        + parameters.control.b_p * com_from_step_midstance ...
                        + parameters.control.b_d * com_vel_midstance ...
                        ;
            elseif strcmp(mode, 'alternating')
                % Equation 4 (eq_footPlacement_ml)
                p_new = p ...
                        + (-1)^n * parameters.control.b_o ...
                        + parameters.control.b_p * com_from_step_midstance ...
                        + parameters.control.b_d * com_vel_midstance ...
                        ;
            else
                error('Variable "mode" must be either ''progressive'' or ''alternating''')
            end
            
            % update discrete variables
            n = n+1;
            step_duration_n = [step_duration_n; parameters.control.T_step]; %#ok<AGROW>
            step_time_n = [step_time_n; t]; %#ok<AGROW>
            foot_placement_n = [foot_placement_n, p_new - p]; %#ok<AGROW>
            this_step_distance = (p_new - p);
            step_distance_n = [step_distance_n; this_step_distance]; %#ok<AGROW>
            com_from_cop_at_step_end_n = [com_from_cop_at_step_end_n, x - p]; %#ok<AGROW>
            com_vel_at_step_end_n = [com_vel_at_step_end_n, v]; %#ok<AGROW>
            eta_step_n = [eta_step_n, eta_step]; %#ok<AGROW>
            p = p_new;
            
            % store indices
            step_indices_discrete = [step_indices_discrete n]; %#ok<AGROW>
            step_indices_continuous = [step_indices_continuous i_time]; %#ok<AGROW>
        end
        
        % Euler step
        a = omega^2 * (x - p); % Equation 1 (eq_biomechanics)
        x = x + dt*v + dt^2*a;
        v = v + dt*a;
                    
        % output
        p_timeseries(i_time) = p;
        x_timeseries(i_time) = x;
        v_timeseries(i_time) = v;
        a_timeseries(i_time) = a;
        u_timeseries(i_time) = 0;
        if strcmp(mode, 'progressive')
            q_timeseries(i_time) = x - p;
        elseif strcmp(mode, 'alternating')
            q_timeseries(i_time) = (x - p)*(-1)^n;
        end
        xcom_timeseries(i_time) = x + 1/omega * v;
        step_number_timeseries(i_time) = n;
        
        if abs(p_timeseries(i_time)) > 50
            break
        end
    end    
    
% ----------------------------------------------------------------------------------------------------------------------
% package results
% ----------------------------------------------------------------------------------------------------------------------
    
    % calculate averages after relaxation time
    com_from_stance_ankle_data_mean_free_left = (com_from_cop_at_midstance_n - mean(com_from_cop_at_midstance_n))';
    com_vel_heelstrike_data_mean_free_left = (com_vel_at_midstance_n - mean(com_vel_at_midstance_n))';
    com_from_stance_ankle_data_mean_free_right = (com_from_cop_at_midstance_n - mean(com_from_cop_at_midstance_n))';
    com_vel_heelstrike_data_mean_free_right = (com_vel_at_midstance_n - mean(com_vel_at_midstance_n))';
    
    foot_placement_x_data_mean_free_left = (foot_placement_n - mean(foot_placement_n))';
    foot_placement_x_data_mean_free_right = (foot_placement_n - mean(foot_placement_n))';    
    
    % package outputs
    results = struct;
    results.continuous = struct;
    results.continuous.time = time;
    results.continuous.p = p_timeseries;
    results.continuous.x = x_timeseries;
    results.continuous.v = v_timeseries;
    results.continuous.a = a_timeseries;
    results.continuous.u = u_timeseries;
    results.continuous.q = q_timeseries;
    results.continuous.xcom = xcom_timeseries;
    results.continuous.step_number = step_number_timeseries;
    results.continuous.metronome = metronome_timeseries;
    results.continuous.metronome_phase = metronome_phase;
    results.continuous.metronome_reset_indices = metronome_reset_indices;
    
    results.discrete = struct;
    results.discrete.step_duration = step_duration_n;
    results.discrete.step_time = step_time_n;
    results.discrete.foot_placement = foot_placement_n;
    results.discrete.step_distance = step_distance_n;
    results.discrete.com_from_stance_ankle_data_mean_free_left = com_from_stance_ankle_data_mean_free_left;
    results.discrete.com_vel_heelstrike_data_mean_free_left = com_vel_heelstrike_data_mean_free_left;
    results.discrete.foot_placement_x_data_mean_free_left = foot_placement_x_data_mean_free_left;
    results.discrete.com_from_stance_ankle_data_mean_free_right = com_from_stance_ankle_data_mean_free_right;
    results.discrete.com_vel_heelstrike_data_mean_free_right = com_vel_heelstrike_data_mean_free_right;
    results.discrete.foot_placement_x_data_mean_free_right = foot_placement_x_data_mean_free_right;
    results.discrete.x_n = results.continuous.x(midstance_indices);
    results.discrete.q_n = results.continuous.q(midstance_indices);
    results.discrete.v_n = results.continuous.v(midstance_indices);
    
    results.indices = struct;
    results.indices.step_indices_continuous = step_indices_continuous;
    results.indices.left_step_indices_continuous = left_step_indices_continuous;
    results.indices.right_step_indices_continuous = right_step_indices_continuous;
    results.indices.step_indices_discrete = step_indices_discrete;
    results.indices.left_step_indices_discrete = left_step_indices_discrete;
    results.indices.right_step_indices_discrete = right_step_indices_discrete;
    
    results.model_parameters.omega = omega;
    
    control = parameters.control;
end
