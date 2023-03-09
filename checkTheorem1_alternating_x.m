% define free variables
syms p_n
syms x_n
syms v_n
syms b_o
syms b_p
syms b_d
syms omega
syms s
syms c
syms n
evalin(symengine,'assume(n,Type::Integer)')

% intermediate variables and identities
p_np1 = p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n;
x_np1 = p_n + (-1)^n*b_o*(1-c) + (c^2 + s^2 + b_p*(1-c))*(x_n-p_n) + (b_d*(1-c) + 2*c*s/omega)*v_n; % equation XX from Lemma 1 (eq_transition_ap)
v_np1 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*(x_n-p_n) + (c^2 + s^2 - b_d*s*omega)*v_n;       % equation XX from Lemma 1 (eq_transition_ap)

% define hyperbolic identities to be used
hyp_id_1 = c^2 + s^2 == 2*s^2 + 1;
hyp_id_2 = c^2 - s^2 == 1;
hyp_id_3 = c^2 + s^2 == 2*c^2 - 1;

% ----------------------------------------------------------------------------------------------------------------------
% alternating
% ----------------------------------------------------------------------------------------------------------------------

% set up equations
eq_1 = (x_n - p_n) == -(x_np1 - p_np1);
eq_2 = v_n == -v_np1;

solution = solve([eq_1, eq_2], [x_n, v_n]);

% simplify x_n
x_n_alt_solution_1 = simplify(solution.x_n);    % solution
x_n_alt_solution_2 = (p_n + 2*c^2*p_n + c^4*p_n + 2*p_n*s^2 + p_n*s^4 + (-1)^n*b_o*c^3 - 2*c^2*p_n*s^2 - b_p*c*p_n + (-1)^n*b_o*c - b_p*c^3*p_n - b_d*omega*p_n*s + b_p*c*p_n*s^2 - b_d*omega*p_n*s^3 - (-1)^n*b_o*c*s^2 + b_d*c^2*omega*p_n*s)/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s + 2*c^2 + b_p*c*s^2 - b_p*c + s^4 - b_d*omega*s^3 + 2*s^2 - b_d*omega*s + 1); % copy explicit solution as basis for simplification
x_n_alt_solution_3 = ((-1)^n*b_o*c*(c^2 - s^2 + 1) + p_n*((c^2 - s^2)^2 + 1 + 2*(c^2 + s^2) - (c^2 - s^2 + 1)*b_p*c - (- (c^2 - s^2) + 1)*b_d*omega*s)) / ((c^2 - s^2)^2 + 2*(c^2 + s^2) + (-(c^2 - s^2) - 1)*c*b_p + (c^2 - s^2 - 1)*b_d*omega*s + 1); % simplify
x_n_alt_solution_4 = ((-1)^n*b_o*c*(1         + 1) + p_n*((1        )^2 + 1 + 2*(2*c^2 - 1) - (1         + 1)*b_p*c - (- (1        ) + 1)*b_d*omega*s)) / ((1        )^2 + 2*(2*c^2 - 1) + (-(1        ) - 1)*c*b_p + (1         - 1)*b_d*omega*s + 1); % plug in hyperbolic identity
x_n_alt_solution_5 = ((-1)^n*b_o + p_n*(2*c - b_p)) / (2*c - b_p); % simplify
x_n_alt_solution_6 = p_n + (-1)^n*b_o / (2*c - b_p); % simplify
% => the solution is q_n = x_n - p_n = (-1)^n*b_o / (2*c - b_p)

disp(simplify(x_n_alt_solution_1 - x_n_alt_solution_2))
disp(simplify(x_n_alt_solution_1 - x_n_alt_solution_3))
disp(simplify(x_n_alt_solution_4 - x_n_alt_solution_5))

% simplify v_n
v_n_alt_solution_1 = simplify(solution.v_n);    % solution
v_n_alt_solution_2 = ((-1)^n*b_o*omega*s*(- c^2 + s^2 + 1))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s + 2*c^2 + b_p*c*s^2 - b_p*c + s^4 - b_d*omega*s^3 + 2*s^2 - b_d*omega*s + 1);    % explicit copy as basis for simplification
v_n_alt_solution_3 = ((-1)^n*b_o*omega*s*(- (c^2 - s^2) + 1))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s + 2*c^2 + b_p*c*s^2 - b_p*c + s^4 - b_d*omega*s^3 + 2*s^2 - b_d*omega*s + 1);    % explicit copy as basis for simplification
v_n_alt_solution_4 = ((-1)^n*b_o*omega*s*(- (1        ) + 1))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s + 2*c^2 + b_p*c*s^2 - b_p*c + s^4 - b_d*omega*s^3 + 2*s^2 - b_d*omega*s + 1);    % use hyperbolic identity
% the numerator is 0, so the solution is v_n == 0

disp(simplify(v_n_alt_solution_1 - v_n_alt_solution_2))
disp(simplify(v_n_alt_solution_1 - v_n_alt_solution_3))











