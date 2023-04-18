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
xt1 = p_n + (x_n - p_n)*c + v_n/omega * s;              % Equation 17
vt1 = (x_n - p_n)*s*omega + v_n*c;                      % Equation 17
p_np1 = p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n;     % Equation 4

% Equation 17 is application of Equation 2
% Equation 18 is plugging in the shorthand definition in Equation 5

disp('This code simplifies equations and then calculates the differences between the old and the new version.')
disp('Differences are printed into the command line and should be 0.')

% q_1 - alternating
% this is the state of x at midstance 1
x_np1_1 = p_np1 + (xt1 - p_np1)*c + vt1*s/omega;    % Equation 2
x_np1_2 = (p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n) + ((p_n + (x_n - p_n)*c + v_n/omega * s) - (p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n))*c + ((x_n - p_n)*s*omega + v_n*c)*s/omega;    % plug in definition of xt1, vt1 and p_np1 from above
x_np1_3 = p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n + p_n*c + (x_n - p_n)*c^2 + v_n/omega*s*c - p_n*c - (-1)^n*b_o*c - b_p*(x_n-p_n)*c - b_d*v_n*c + (x_n - p_n)*s^2 + v_n*c*s/omega;    % expand
x_np1_4 = p_n + (-1)^n*b_o*(1-c) + (c^2 + s^2 + b_p*(1-c))*(x_n-p_n) + (b_d*(1-c) + 2*c*s/omega)*v_n;    % sort
disp(simplify(x_np1_1 - x_np1_2))
disp(simplify(x_np1_1 - x_np1_3))
disp(simplify(x_np1_1 - x_np1_4))

% v_1 - alternating
v_np1_1 = (xt1 - p_np1)*s*omega + vt1*c; % Equation 2
v_np1_2 = ((p_n + (x_n - p_n)*c + v_n/omega * s) - (p_n + (-1)^n*b_o + b_p*(x_n-p_n) + b_d*v_n))*s*omega + ((x_n - p_n)*s*omega + v_n*c)*c; % plug in definition of xt1, vt1 and p_np1 from above
v_np1_3 = (x_n - p_n)*c*s*omega + v_n*s^2 - (-1)^n*b_o*s*omega - b_p*(x_n-p_n)*s*omega - b_d*v_n*s*omega + (x_n - p_n)*s*c*omega + v_n*c^2; % expand
v_np1_4 = - (-1)^n*b_o*s*omega + (2*c - b_p)*s*omega*(x_n-p_n) + (c^2 + s^2 - b_d*s*omega)*v_n; % sort
disp(simplify(v_np1_1 - v_np1_2))
disp(simplify(v_np1_1 - v_np1_3))
disp(simplify(v_np1_1 - v_np1_4))

