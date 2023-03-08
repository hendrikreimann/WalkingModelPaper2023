% define free variables
syms b_o
syms b_p
syms b_d
syms omega
syms s
syms c
syms n
evalin(symengine,'assume(n,Type::Integer)')
syms q_n
syms v_n

% intermediate variables and identities
q_np1_pro = - c*b_o + (c^2 + s^2 - c*b_p)*q_n - (b_d - 2*s/omega)*c*v_n;                                        % equation XX from Lemma 1 (eq_transition_ap)
v_np1_pro = - s*omega*b_o + (2*c - b_p)*s*omega*q_n - (s*omega*b_d - c^2 - s^2)*v_n;                            % equation XX from Lemma 1 (eq_transition_ap)

q_np1_alt = c*b_o + (-1)^(n+1)*(c^2 + s^2 - c*b_p)*q_n - (-1)^(n+1)*(b_d - 2*s/omega)*c*v_n;                    % equation XX from Lemma 1 (eq_transition_ml)
v_np1_alt = -(-1)^n*s*omega*b_o + (-1)^(n)*(2*c - b_p)*s*omega*(-1)^(n)*q_n - (b_d*s*omega - c^2 - s^2)*v_n;    % equation XX from Lemma 1 (eq_transition_ml)

% define hyperbolic identities to be used
hyp_id_1 = c^2 + s^2 == 2*s^2 + 1;
hyp_id_2 = c^2 - s^2 == 1;
hyp_id_3 = c^2 + s^2 == 2*c^2 - 1;

% ----------------------------------------------------------------------------------------------------------------------
% alternating
% ----------------------------------------------------------------------------------------------------------------------

% set up equations
eq_1_alt = q_n == - q_np1_alt;
eq_2_alt = v_n == - v_np1_alt;

% solve
solution_alt = solve([eq_1_alt, eq_2_alt], [q_n, v_n]);




% simplify solution for q_0 using hyperbolic identities
q_0_alt_solution_1 = simplify(solution_alt.q_n);    % solution
q_0_alt_solution_2 = -(b_o*c*(c^2 - s^2 + 1))/(c^2 + s^2 - (-1)^n*c^2 - (-1)^n*c^4 - (-1)^n*s^2 - (-1)^n*s^4 + (-1)^n*b_p*c^3 - b_d*omega*s + 2*(-1)^n*c^2*s^2 + (-1)^n*b_p*c - (-1)^n*b_p*c*s^2 + (-1)^n*b_d*omega*s^3 - (-1)^n*b_d*c^2*omega*s + 1);  % paste explicit solution as basis for simplification
q_0_alt_solution_3 = -(b_o*c*(c^2 - s^2 + 1)) / (- (-1)^n*(c^2 - s^2)^2 + (c^2 + s^2) - (-1)^n*(c^2 + s^2) + (-1)^n*(c^2 - s^2 + 1)*c*b_p - (1 + (-1)^n*(c^2 - s^2))*b_d*omega*s + 1);  % simplify
q_0_alt_solution_4 = -(b_o*c*(1         + 1)) / (- (-1)^n*(1        )^2 + (2*s^2 + 1) - (-1)^n*(2*s^2 + 1) + (-1)^n*(1         + 1)*c*b_p - (1 + (-1)^n*(1        ))*b_d*omega*s + 1);  % plug in hyperbolic identities
q_0_alt_solution_5 = -(2*b_o*c) / ( - (-1)^n*2 + 2*s^2 - (-1)^n*2*s^2 + (-1)^n*2*c*b_p - (1 + (-1)^n)*b_d*omega*s + 2);  % simplify

disp(simplify(q_0_alt_solution_1 - q_0_alt_solution_2))
disp(simplify(q_0_alt_solution_1 - q_0_alt_solution_3))
disp(simplify(q_0_alt_solution_4 - q_0_alt_solution_5))

pretty(q_0_alt_solution_5)



return
q_0_alt_solution_3 = -(b_o*c)/(c^2 - s^2 - b_p*c + b_d*omega*s + 1);  % simplify




disp(simplify(q_0_alt_solution_1 - q_0_alt_solution_3))

% return
% ----------------------------------------------------------------------------------------------------------------------
% progressive
% ----------------------------------------------------------------------------------------------------------------------

% set up equations
eq_1_pro = q_n == q_np1_pro;
eq_2_pro = v_n == v_np1_pro;

% solve
solution_pro = solve([eq_1_pro, eq_2_pro], [q_n, v_n]);

% simplify solution for q_0 using hyperbolic identities
q_n_pro_solution_1 = simplify(solution_pro.q_n);    % solution
q_n_pro_solution_2 = -(b_o*c*( - c^2 + s^2 + 1))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s - 2*c^2 + b_p*c*s^2 + b_p*c + s^4 - b_d*omega*s^3 - 2*s^2 + b_d*omega*s + 1);    % paste explicit solution as basis for simplification
q_n_pro_solution_3 = -(b_o*c*(-(c^2 - s^2 - 1)))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s - 2*c^2 + b_p*c*s^2 + b_p*c + s^4 - b_d*omega*s^3 - 2*s^2 + b_d*omega*s + 1);    % simplify
q_n_pro_solution_4 = -(b_o*c*(-(1         - 1)))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s - 2*c^2 + b_p*c*s^2 + b_p*c + s^4 - b_d*omega*s^3 - 2*s^2 + b_d*omega*s + 1);    % apply hyporbolic identity hyp_id_2
% the numerator is 0, so the solution is q_0 == 0

% verify transformations
disp(simplify(q_n_pro_solution_1 - q_n_pro_solution_2))
disp(simplify(q_n_pro_solution_1 - q_n_pro_solution_3))

% simplify solution for v_0 using hyperbolic identities
v_n_pro_solution_1 = simplify(solution_pro.v_n);    % solution
v_n_pro_solution_2 = -(b_o*omega*s*(c^2 - s^2 + 1))/(c^4 - b_p*c^3 - 2*c^2*s^2 + b_d*omega*c^2*s - 2*c^2 + b_p*c*s^2 + b_p*c + s^4 - b_d*omega*s^3 - 2*s^2 + b_d*omega*s + 1); % paste explicit solution as basis for simplification
v_n_pro_solution_3 = -(b_o*omega*s*(c^2 - s^2 + 1))/((c^2 - s^2)^2 - 2*(c^2 + s^2) - b_p*c*(c^2 - s^2 - 1) + b_d*omega*s*(c^2 - s^2 + 1) + 1); % simplify
v_n_pro_solution_4 = -(b_o*omega*s*(1         + 1))/((1        )^2 - 2*(2*s^2 + 1) - b_p*c*(1         - 1) + b_d*omega*s*(1         + 1) + 1); % apply hyporbolic identity hyp_id_2
v_n_pro_solution_5 = -(b_o*omega*s*(1         + 1))/((1        )^2 - 2*(2*s^2 + 1) - b_p*c*(1         - 1) + b_d*omega*s*(1         + 1) + 1); % simplify
v_n_pro_solution_6 = (b_o*omega)/(2*s - b_d*omega); % simplify
% this is the solution in Theorem 1

% verify transformations
disp(simplify(v_n_pro_solution_1 - v_n_pro_solution_2))
disp(simplify(v_n_pro_solution_1 - v_n_pro_solution_3))
disp(simplify(v_n_pro_solution_4 - v_n_pro_solution_5))
disp(simplify(v_n_pro_solution_4 - v_n_pro_solution_6))







