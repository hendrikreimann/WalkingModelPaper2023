% Script to check Equation 2 (eq_analyticSystemSolution) of the paper.

disp('Solving walking system for one very long stance duration with random initial contitions');
disp('This table compares the numeric with the analytical solution from Equation 2 in the paper:')

% set parameters
dt = 1e-6;                  % very low Euler step, to maximize accuracy
stepping_frequency = 0.1;   % this is set to a low number to get more time before the system takes a step (which is a different equation)

% set parameters
parameters.biomechanics = struct;
T_total = 1/stepping_frequency * 0.25;
T_step = 1/(2*stepping_frequency);
parameters.biomechanics.z_c = 0.814;            % repeats the values in the simulation
parameters.biomechanics.g = 9.81;
omega = sqrt(parameters.biomechanics.g/parameters.biomechanics.z_c);

% define reference configuration
p = rand;
x_0 = rand;
v_0 = rand;

% define control struct
parameters.control = struct;
parameters.control.T_step = T_step;
parameters.control.b_o = 0;
parameters.control.b_p = 0;
parameters.control.b_d = 0;

% calculate numeric solutions
results = simulateSlipWalker(T_total, parameters, x_0, v_0, p, dt);
x_end_numeric = results.continuous.x(end);
v_end_numeric = results.continuous.v(end);
p_end_numeric = results.continuous.p(end);

% calculate analytic solutions for the same time
time = results.continuous.time;
x_analytic = p + (x_0 - p)*cosh(omega * time) + v_0/omega*sinh(omega * time);
v_analytic = (x_0 - p)*sinh(omega * time)*omega + v_0*cosh(omega * time);

% collect results and print to command line
x_end_analytic = x_analytic(end);
v_end_analytic = v_analytic(end);
x_end_error = abs(x_end_numeric-x_end_analytic);
v_end_error = abs(v_end_numeric-v_end_analytic);
results_table = table ...
  ( ...
    [x_end_numeric; v_end_numeric], [x_end_analytic; v_end_analytic], [x_end_error; v_end_error], ...
    'VariableNames', {'numeric', 'analytic', 'error'}, ...
    'RowNames', {'x', 'v'} ...
  );
disp(results_table);

% plot results
figure; tiledlayout(2, 1, 'TileSpacing', 'Compact', 'padding', 'compact');

nexttile; hold on
plot(results.continuous.time, results.continuous.x, 'linewidth', 4, 'DisplayName', 'CoM numeric');
plot(time, x_analytic, '-', 'linewidth', 2, 'DisplayName', 'CoM analytic');
legend('Location', 'best');
xlabel('time (s)', 'fontsize', 14)
ylabel('lateral position', 'fontsize', 14)

nexttile; hold on
plot(results.continuous.time, results.continuous.v, 'linewidth', 4, 'DisplayName', 'CoM numeric');
plot(time, v_analytic, '-', 'linewidth', 2, 'DisplayName', 'CoM analytic');
legend('Location', 'best');
xlabel('time (s)', 'fontsize', 14)
ylabel('velocity', 'fontsize', 14)
    




