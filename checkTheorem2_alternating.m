% define free variables
syms b_o
syms b_p
syms b_d
syms omega
syms s
syms c
syms n
evalin(symengine,'assume(n,Type::Integer)')
syms p_n
syms x_n
syms v_n

% define dependent variables
q_n = x_n - p_n;                                                                                        % definition
q_n_ref = (-1)^n*b_o/(2*c - b_p);                                                                       % from Theorem 1
v_n_ref = 0;                                                                                            % from Theorem 1
x_np1 = p_n + (-1)^n*b_o*(1-c) + (c^2 + s^2 + b_p*(1-c))*(x_n-p_n) + (b_d*(1-c) + 2*c*s/omega)*v_n;     % from Lemma 1
p_np1 = p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n;                                                     % Equation 4 (eq_footPlacement_ml)
v_np1 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*(x_n-p_n) + (c^2 + s^2 - b_d*s*omega)*v_n;           % from Lemma 1
q_np1 = x_np1 - p_np1;                                                                                  % definition
q_np1_ref = (-1)^(n+1)*b_o / (2*c - b_p);                                                               % from Theorem 1

delta_q_n = q_n - q_n_ref;                                                                              % definition
delta_v_n = v_n - v_n_ref;                                                                              % definition

% define hyperbolic identities to be used
hyp_id_3 = c^2 + s^2 == 2*c^2 - 1;

% calculate new delta_q
delta_q_np1_1 = q_np1 - q_np1_ref;                                                % definition
delta_q_np1_2 = x_np1 - p_np1 - (-1)^(n+1)*b_o / (2*c - b_p);                       % plug in equations from above
delta_q_np1_3 = p_n + (-1)^n*b_o*(1-c) + (c^2 + s^2 + b_p*(1-c))*(x_n-p_n) + (b_d*(1-c) + 2*c*s/omega)*v_n - (p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n) - (-1)^(n+1)*b_o / (2*c - b_p);                       % plug in more equations from above
delta_q_np1_4 = - (-1)^n*b_o*c + (c^2 + s^2 - b_p*c)*q_n + (2*s/omega - b_d)*c*v_n - (-1)^(n+1)*b_o / (2*c - b_p);                       % expand and simplify
delta_q_np1_5 = - (-1)^n*b_o*c + (c^2 + s^2 - b_p*c)*(delta_q_n + q_n_ref) + (2*s/omega - b_d)*c*(delta_v_n + v_n_ref) - (-1)^(n+1)*b_o / (2*c - b_p);                       % plug in q_n and v_n
delta_q_np1_6 = - (-1)^n*b_o*c + (c^2 + s^2 - b_p*c)*(delta_q_n + (-1)^n*b_o / (2*c - b_p)) + (2*s/omega - b_d)*c*delta_v_n - (-1)^(n+1)*b_o / (2*c - b_p);                       % plug in q_n_ref and v_n_ref
delta_q_np1_7 = - (-1)^n*(b_o*c*(2*c - b_p) - b_o - (c^2 + s^2 - b_p*c)*b_o)/(2*c - b_p) + (c^2 + s^2 - b_p*c)*delta_q_n + (2*s/omega - b_d)*c*delta_v_n;                       % expand and sort
delta_q_np1_8 = - (-1)^n*(b_o*c*(2*c - b_p) - b_o - (2*c^2 - 1 - b_p*c)*b_o)/(2*c - b_p) + (c^2 + s^2 - b_p*c)*delta_q_n + (2*s/omega - b_d)*c*delta_v_n;                       % use hyp_id_3
delta_q_np1_9 = (c^2 + s^2 - c*b_p)*delta_q_n - (b_d - 2*s/omega)*c*delta_v_n;                       % simplify

% verify transformations
disp(simplify(delta_q_np1_1 - delta_q_np1_2))
disp(simplify(delta_q_np1_1 - delta_q_np1_3))
disp(simplify(delta_q_np1_1 - delta_q_np1_4))
disp(simplify(delta_q_np1_1 - delta_q_np1_5))
disp(simplify(delta_q_np1_1 - delta_q_np1_6))
disp(simplify(delta_q_np1_1 - delta_q_np1_7))
disp(simplify(delta_q_np1_8 - delta_q_np1_9))

% calculate new delta_v
delta_v_np1_1 = v_np1 - v_n_ref;                                                    % definition
delta_v_np1_2 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*q_n + (c^2 + s^2 - b_d*s*omega)*v_n;   % plug in equations from above
delta_v_np1_3 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*(delta_q_n + q_n_ref) + (c^2 + s^2 - b_d*s*omega)*(delta_v_n + v_n_ref);   % plug in q_n and v_n
delta_v_np1_4 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*(delta_q_n + (-1)^n*b_o/(2*c - b_p)) + (c^2 + s^2 - b_d*s*omega)*delta_v_n;   % plug in q_n_ref and v_n_ref
delta_v_np1_5 = - (-1)^n*(b_o*s*omega - (2*c - b_p)*s*omega*b_o/(2*c - b_p)) + (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - b_d*s*omega)*delta_v_n;   % expand
delta_v_np1_6 = (2*c - b_p)*s*omega*delta_q_n + (c^2 + s^2 - s*omega*b_d)*delta_v_n; % simplify


% verify transformations
disp(simplify(delta_v_np1_1 - delta_v_np1_2))
disp(simplify(delta_v_np1_1 - delta_v_np1_3))
disp(simplify(delta_v_np1_1 - delta_v_np1_4))
disp(simplify(delta_v_np1_1 - delta_v_np1_5))
disp(simplify(delta_v_np1_1 - delta_v_np1_6))













