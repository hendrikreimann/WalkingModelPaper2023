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

% define dependent variables
q_n_ref = 0;                                                                        % from Theorem 1
v_n_ref = b_o/(2*s/omega - b_d);                                                    % from Theorem 1
q_np1 = - c*b_o + (c^2 + s^2 - c*b_p)*q_n - (b_d - 2*s/omega)*c*v_n;                % from Lemma 1
v_np1 = - s*omega*b_o + (2*c - b_p)*s*omega*q_n + (c^2 + s^2 - s*omega*b_d)*v_n;    % from Lemma 1
delta_q_n = q_n - q_n_ref;                                                          % definition
delta_v_n = v_n - v_n_ref;                                                          % definition

% define hyperbolic identities to be used
hyp_id_1 = c^2 + s^2 == 2*s^2 + 1;

% calculate new delta_q
delta_q_np1_1 = q_np1 - q_n_ref;                                                % definition
delta_q_np1_2 = - c*b_o + (c^2 + s^2 - c*b_p)*q_n - (b_d - 2*s/omega)*c*v_n;    % plug in equations from above
delta_q_np1_3 = - c*b_o + (c^2 + s^2 - c*b_p)*(delta_q_n + q_n_ref) - (b_d - 2*s/omega)*c*(delta_v_n + v_n_ref);    % plug in q_n and v_n
delta_q_np1_4 = - c*b_o + (c^2 + s^2 - c*b_p)*delta_q_n - (b_d - 2*s/omega)*c*(delta_v_n + b_o/(2*s/omega - b_d));    % plug in q_n_ref and v_n_ref
delta_q_np1_5 = - c*b_o + (c^2 + s^2 - c*b_p)*delta_q_n - (b_d - 2*s/omega)*c*delta_v_n  - (b_d - 2*s/omega)*c*b_o/(2*s/omega - b_d);    % expand
delta_q_np1_6 = (c^2 + s^2 - c*b_p)*delta_q_n - (b_d - 2*s/omega)*c*delta_v_n;    % simplify

% verify transformations
disp(simplify(delta_q_np1_1 - delta_q_np1_2))
disp(simplify(delta_q_np1_1 - delta_q_np1_3))
disp(simplify(delta_q_np1_1 - delta_q_np1_4))
disp(simplify(delta_q_np1_1 - delta_q_np1_5))
disp(simplify(delta_q_np1_1 - delta_q_np1_6))

% calculate new delta_v
delta_v_np1_1 = v_np1 - v_n_ref;                                                    % definition
delta_v_np1_2 = - s*omega*b_o + (2*c - b_p)*s*omega*q_n + (c^2 + s^2 - s*omega*b_d)*v_n - b_o/(2*s/omega - b_d); % plug in equations from above
delta_v_np1_3 = - s*omega*b_o + (2*c - b_p)*s*omega*(delta_q_n + q_n_ref) + (c^2 + s^2 - s*omega*b_d)*(delta_v_n + v_n_ref) - b_o/(2*s/omega - b_d); % plug in q_n and v_n
delta_v_np1_4 = - s*omega*b_o + (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - s*omega*b_d)*(delta_v_n + b_o/(2*s/omega - b_d)) - b_o/(2*s/omega - b_d); % plug in q_n_ref and v_n_ref
delta_v_np1_5 = - s*omega*b_o + (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - s*omega*b_d)*delta_v_n + (c^2 + s^2 - s*omega*b_d)*b_o/(2*s/omega - b_d) - b_o/(2*s/omega - b_d); % expand
delta_v_np1_6 = - s*omega*b_o + (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - s*omega*b_d)*delta_v_n + (2*s^2 + 1 - s*omega*b_d)*b_o/(2*s/omega - b_d) - b_o/(2*s/omega - b_d); % plug in hyperbolic identity 
delta_v_np1_7 = (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - s*omega*b_d)*delta_v_n; % simplify

% verify transformations
disp(simplify(delta_v_np1_1 - delta_v_np1_2))
disp(simplify(delta_v_np1_1 - delta_v_np1_3))
disp(simplify(delta_v_np1_1 - delta_v_np1_4))
disp(simplify(delta_v_np1_1 - delta_v_np1_5))
disp(simplify(delta_v_np1_6 - delta_v_np1_7))

