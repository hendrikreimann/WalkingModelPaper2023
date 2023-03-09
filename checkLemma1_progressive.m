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

% define intermediate variables
xt1 = p_n + (x_n - p_n)*c + v_n/omega * s;              % Equation 24-26 (eq_lemma1_1 - 3)
vt1 = (x_n - p_n)*s*omega + v_n*c;                      % Equation 24-26 (eq_lemma1_1 - 3)
p_np1 = p_n + b_o + b_p*(x_n-p_n) + b_d*v_n;            % Equation 3 (eq_footPlacement_ap)
x_np1 = p_np1 + (xt1 - p_np1)*c + vt1*s/omega;          % equation 2 (eq_analyticSystemSolution)
q_n = x_n - p_n;                                        % equation 5 (eq_def_q_prg)

% Equation 24 (eq_lemma1_1) is application of Equation 2 (eq_analyticSystemSolution)
% Equation 25 (eq_lemma1_2) is plugging in the shorthand definition in Equation 5 (eq_hyperbolicConstants)
% Equation 26 (eq_lemma1_3) is plugging in the definition of q in Equation 6 (eq_def_q_prg)

disp('This code simplifies equations and then calculates the differences between the old and the new version.')
disp('Differences are printed into the command line and should be 0.')

% q_1 - progressive
% this is the state of x at midstance n+1
q_np1_1_pro = x_np1 - p_np1;
q_np1_2_pro = p_np1 + (xt1-p_np1)*c + vt1*s/omega - p_np1; % plug in equation for x_1_pro from above
q_np1_3_pro = ((p_n + (x_n - p_n)*c + v_n/omega * s)-(p_n + b_o + b_p*(x_n-p_n) + b_d*v_n))*c + ((x_n - p_n)*s*omega + v_n*c)*s/omega; % plug in definition of xt1, vt1 and p_1_pro from above
q_np1_4_pro = c^2*q_n + c*s/omega*v_n - c*b_o - c*b_p*q_n - c*b_d*v_n + s^2*q_n + c*s/omega*v_n; % expand and plug in equation for q_0_pro from above
q_np1_5_pro = - c*b_o + (c^2 + s^2 - c*b_p)*q_n - (b_d - 2*s/omega)*c*v_n; % sort and simplify
disp(simplify(q_np1_1_pro - q_np1_2_pro))
disp(simplify(q_np1_1_pro - q_np1_3_pro))
disp(simplify(q_np1_1_pro - q_np1_4_pro))
disp(simplify(q_np1_1_pro - q_np1_5_pro))

% v_1 - progressive
v_np1_1_pro = (xt1 - p_np1)*s*omega + vt1*c; % apply equation for analytical solution
v_np1_2_pro = ((p_n + (x_n - p_n)*c + v_n/omega * s) - (p_n + b_o + b_p*(x_n-p_n) + b_d*v_n))*s*omega + ((x_n - p_n)*s*omega + v_n*c)*c; % plug in definition of xt1, vt1 and p_1_pro from above
v_np1_3_pro = q_n*c*s*omega + v_n*s^2 - b_o*s*omega - b_p*q_n*s*omega - b_d*v_n*s*omega + q_n*s*omega*c + v_n*c^2; % expand and plug in equation for q_0_pro from above
v_np1_4_pro = - s*omega*b_o + (2*c - b_p)*s*omega*q_n + (c^2 + s^2 - s*omega*b_d)*v_n; % sort and simplify
disp(simplify(v_np1_1_pro - v_np1_2_pro))
disp(simplify(v_np1_1_pro - v_np1_3_pro))
disp(simplify(v_np1_1_pro - v_np1_4_pro))

































