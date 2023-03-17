% determining preferential cadence for speed of 1.0, 1.2 and 1.4 m/s according to Inman formula
% step_length = stridelength_cadence_ratio * body_height * cadence * 0.5
% with stridelength_cadence_ratio = 0.008, cadence in steps per minute

% Inman, V.T., H.J. Ralston, F. Todd, and J.C. Lieberman. Human Walking. Williams & Wilkins, 1981. 
% Wang, Yang, and Manoj Srinivasan. “Stepping in the Direction of the Fall: The next Foot Placement Can Be Predicted 
%   from Current Upper Body State in Steady-State Walking.” Biology Letters 10 (2014): 20140405.

% variable
speed = 1.2; % in meters / seconds -- change this
body_height = 1.74; % average from Wang and Srinivasan, 2014

% parameter (from Inman)
stridelength_cadence_ratio = 0.008; % normalized to body_height, in seconds
minutes_to_seconds = 60;

%% transform equations
% speed = step_length / step_time;    % definition of speed, in meters / seconds
% speed = step_length * cadence_in_steps_per_second;    % definition of speed, in meters / seconds
% step_length = speed * cadence_in_steps_per_second^(-1);    % step length, in meters
% step_length = speed * (cadence * 60^(-1))^(-1);    % step length, in meters
% cadence_in_steps_per_second = cadence * 60^(-1);

%% we get the stridelength_cadence_ratio from Inman. 
% (2 * step_length/body_height) / cadence = stridelength_cadence_ratio

%% now solve it for cadence and remove step_length = 
% cadence = (2 * step_length/body_height) / stridelength_cadence_ratio; % in steps / minute
% cadence = (2 * step_length/body_height) / stridelength_cadence_ratio; % in steps / minute
% cadence = 2 * step_length * body_height^(-1) * stridelength_cadence_ratio^(-1); % in steps / minute
% cadence = 2 * speed * (cadence * 60^(-1))^(-1) * body_height^(-1) * stridelength_cadence_ratio^(-1); % in steps / minute
% cadence^2 = 2 * speed * body_height^(-1) * 60 * stridelength_cadence_ratio^(-1); % in (steps / minute)^2
% cadence = sqrt(2 * speed * body_height^(-1) * 60 * stridelength_cadence_ratio^(-1)); % in (steps / minute)

% set variables
cadence = sqrt(2 * speed * body_height^(-1) * 60 * stridelength_cadence_ratio^(-1)) % in (steps / minute)
