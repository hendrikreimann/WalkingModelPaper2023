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

% define hyperbolic identities to be used
hyp_id_1 = c^2 + s^2 == 2*s^2 + 1;
hyp_id_2 = c^2 - s^2 == 1;
hyp_id_3 = c^2 + s^2 == 2*c^2 - 1;

% define dependent variables
A = [(c^2 + s^2 - c*b_p), (2*s/omega - b_d)*c; (2*c - b_p)*s*omega, (c^2 + s^2 - s*omega*b_d)];

% ----------------------------------------------------------------------------------------------------------------------
% determinant
% ----------------------------------------------------------------------------------------------------------------------
det_A_1 = A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1); % definition of the determinant
det_A_2 = (c^2 + s^2 - c*b_p)*(c^2 + s^2 - s*omega*b_d) - (2*s/omega - b_d)*c*(2*c - b_p)*s*omega; % plug in definition from above
det_A_3 = c^4 + 2*c^2*s^2 - b_p*c^3 + s^4 - c*s^2*b_p - c^2*s*omega*b_d - s^3*omega*b_d + c*s*omega*b_p*b_d - 4*c^2*s^2 + 2*c^2*s*omega*b_d + 2*c*s^2*b_p - c*s*omega*b_d*b_p; % expand
det_A_4 = (c^2 - s^2)^2 - (c^2 - s^2)*c*b_p + (c^2 - s^2)*s*omega*b_d; % sort
det_A_5 = 1             - 1          *c*b_p + 1*          s*omega*b_d; % apply hyperbolic identity

% verify transformations
% disp(simplify(det_A_1 - det_A_2))
% disp(simplify(det_A_1 - det_A_3))
% disp(simplify(det_A_1 - det_A_4))

% ----------------------------------------------------------------------------------------------------------------------
% trace
% ----------------------------------------------------------------------------------------------------------------------

trace_A_1 = A(1, 1) + A(2, 2);
trace_A_2 = (c^2 + s^2 - c*b_p) + (c^2 + s^2 - s*omega*b_d);    % plug in definition from above
trace_A_3 = c^2 + s^2 - c*b_p + c^2 + s^2 - s*omega*b_d;        % expand
trace_A_4 = 2*(c^2 + s^2) - (c*b_p + s*omega*b_d);              % sort
trace_A   = 2*c^2 + 2*s^2 - c*b_p - s*omega*b_d;

% verify transformations
% disp(simplify(trace_A_1 - trace_A_2))
% disp(simplify(trace_A_1 - trace_A_3))
% disp(simplify(trace_A_1 - trace_A_4))
% disp(simplify(trace_A_1 - trace_A))

% ----------------------------------------------------------------------------------------------------------------------
% Discriminant
% ----------------------------------------------------------------------------------------------------------------------

discriminant_A_1 = trace_A_4^2 - 4*det_A_5;                                                     % definition of the discriminant
discriminant_A_2 = (2*(c^2 + s^2) - (c*b_p + s*omega*b_d))^2 - 4*(1 - c*b_p + s*omega*b_d);     % plug in equations from above
discriminant_A_3 = (2*c^2 + 2*s^2 - c*b_p - s*omega*b_d)*(2*c^2 + 2*s^2 - c*b_p - s*omega*b_d) - 4 + 4*c*b_p - 4*s*omega*b_d;     % expand
discriminant_A_4 = 4*c^4 + 4*s^2*c^2 - 2*c^3*b_p - 2*c^2*s*omega*b_d ...
                   + 4*c^2*s^2 + 4*s^4 - 2*c*s^2*b_p - 2*s^3*omega*b_d ...
                   - 2*c^3*b_p - 2*c*s^2*b_p + c^2*b_p^2 + c*s*omega*b_p*b_d ...
                   - 2*c^2*s*omega*b_d - 2*s^3*omega*b_d + c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 ...
                   - 4 + 4*c*b_p - 4*s*omega*b_d;     % expand
discriminant_A_5 = 4*(c^2 + s^2)^2 - 4*(c^2 + s^2 - 1)*c*b_p - 4*(c^2 + s^2 + 1)*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;     % sort
discriminant_A_6 = 4*(c^2 + s^2)^2 - 4*(2*s^2 + 1 - 1)*c*b_p - 4*(2*c^2 - 1 + 1)*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;     % apply hyperbolic identities
discriminant_A_7   = 4*(c^2 + s^2)^2 - 8*s^2*c*b_p - 8*c^2*s*omega*b_d + (c*b_p + s*omega*b_d)^2 - 4;     % simplify
discriminant_A     = 4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p - 8*c^2*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;

% verify transformations
% disp(simplify(discriminant_A_1 - discriminant_A_2))
% disp(simplify(discriminant_A_1 - discriminant_A_3))
% disp(simplify(discriminant_A_1 - discriminant_A_4))
% disp(simplify(discriminant_A_1 - discriminant_A_5))
% disp(simplify(discriminant_A_6 - discriminant_A_7))
% disp(simplify(discriminant_A_6 - discriminant_A))

% ----------------------------------------------------------------------------------------------------------------------
% Eigenvalue - discriminant positive
% ----------------------------------------------------------------------------------------------------------------------

% define Eigenvalue equation
eq_eigenvalues_dpos_1 = trace_A^4 - 2*trace_A^2*discriminant_A - 8*trace_A^2  + (discriminant_A-4)^2;
% expand
eq_eigenvalues_dpos_2 = expand(eq_eigenvalues_dpos_1); 
% copy the result as a basis for futher simplification
eq_eigenvalues_dpos_3 = 16*b_d^2*c^4*omega^2*s^2 - 32*b_d^2*c^2*omega^2*s^4 + 16*b_d^2*omega^2*s^6 - 16*b_d^2*omega^2*s^2 - 32*b_d*b_p*c^5*omega*s + 64*b_d*b_p*c^3*omega*s^3 - 32*b_d*b_p*c*omega*s^5 - 32*b_d*b_p*c*omega*s + 128*b_d*c^2*omega*s + 16*b_p^2*c^6 - 32*b_p^2*c^4*s^2 + 16*b_p^2*c^2*s^4 - 16*b_p^2*c^2 + 128*b_p*c*s^2 - 64*c^4 - 128*c^2*s^2 - 64*s^4 + 64;
% reorder
eq_eigenvalues_dpos_4 = - 64*(c^2 + s^2)^2 + 128*b_p*c*s^2 + 128*b_d*c^2*omega*s + 16*b_p^2*c^2*((c^2 - s^2)^2 - 1) - 32*b_d*b_p*omega*c*s*((c^2 - s^2)^2 + 1) + 16*b_d^2*s^2*omega^2*((c^2 - s^2)^2 - 1) + 64;
% plug in hyperbolic identity
eq_eigenvalues_dpos_5 = - 64*(2*c^2 - 1)^2 + 128*b_p*c*s^2 + 128*b_d*c^2*omega*s + 16*b_p^2*c^2*((1        )^2 - 1) - 32*b_d*b_p*omega*c*s*((1        )^2 + 1) + 16*b_d^2*s^2*omega^2*((1        )^2 - 1) + 64;
% simplify
eq_eigenvalues_dpos_6 = 64*(- 4*c^4 + 4*c^2 + 2*b_p*c*s^2 + 2*b_d*c^2*omega*s - b_d*b_p*omega*c*s);
% simplify further
eq_eigenvalues_dpos_7 = 64*(- 4*c^2*(c^2 - 1) + 2*b_p*c*s^2 + 2*b_d*c^2*omega*s - b_d*b_p*omega*c*s);
% hyperbolic identity
eq_eigenvalues_dpos_8 = 64*(- 4*c^2*(s^2    ) + 2*b_p*c*s^2 + 2*b_d*c^2*omega*s - b_d*b_p*omega*c*s);
% simplify
eq_eigenvalues_dpos_9 = 64*c*s*(2*s - b_d*omega)*(b_p - 2*c);

% verify transformations
% disp(simplify(eq_eigenvalues_dpos_1 - eq_eigenvalues_dpos_2))
% disp(simplify(eq_eigenvalues_dpos_1 - eq_eigenvalues_dpos_3))
% disp(simplify(eq_eigenvalues_dpos_1 - eq_eigenvalues_dpos_4))
% disp(simplify(eq_eigenvalues_dpos_5 - eq_eigenvalues_dpos_6))
% disp(simplify(eq_eigenvalues_dpos_6 - eq_eigenvalues_dpos_7))
% disp(simplify(eq_eigenvalues_dpos_8 - eq_eigenvalues_dpos_9))

% ----------------------------------------------------------------------------------------------------------------------
% Eigenvalue - discriminant negative
% ----------------------------------------------------------------------------------------------------------------------

% define Eigenvalue equation
eq_eigenvalues_dneg_1 = 1/4*trace_A^2 + 1/4*(-discriminant_A) - 1;
% plug in definitions
eq_eigenvalues_dneg_2 = 1/4*(2*c^2 + 2*s^2 - c*b_p - s*omega*b_d)^2 - 1/4*(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p - 8*c^2*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4) - 1;
% expand
eq_eigenvalues_dneg_3 = 1/4*(2*c^2*2*c^2 + 2*s^2*2*c^2 - c*b_p*2*c^2 - s*omega*b_d*2*c^2 + 2*c^2*2*s^2 + 2*s^2*2*s^2 - c*b_p*2*s^2 - s*omega*b_d*2*s^2 - 2*c^2*c*b_p - 2*s^2*c*b_p + c*b_p*c*b_p + s*omega*b_d*c*b_p - 2*c^2*s*omega*b_d - 2*s^2*s*omega*b_d + c*b_p*s*omega*b_d + s*omega*b_d*s*omega*b_d - 4*c^4 - 4*2*c^2*s^2 - 4*s^4 + 8*s^2*c*b_p + 8*c^2*s*omega*b_d - c^2*b_p^2 - 2*c*s*omega*b_p*b_d - s^2*omega^2*b_d^2 + 4 - 4);
% simplify
eq_eigenvalues_dneg_4 = - (c^2 - s^2)*c*b_p + (c^2 - s^2)*s*omega*b_d;
% use hyperbolic identity
eq_eigenvalues_dneg_5 = - (1        )*c*b_p + (1        )*s*omega*b_d;
% simplify
eq_eigenvalues_dneg_6 = s*omega*b_d - c*b_p;

% verify transformations
% disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_2))
% disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_3))
% disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_4))


% ----------------------------------------------------------------------------------------------------------------------
% calculate Eigenvalues for representative points
% ----------------------------------------------------------------------------------------------------------------------

% define origin and directions as basis for the seven candidate points
P_origin = [2*c; 2*s/omega];
d_p = 2*s^2/c - 2*c;
d_v = 2*c^2/(s*omega) - 2*s/omega;
u_p = [d_p; 0];
u_v = [0; d_v];

% point 1
P_1 = P_origin + 1/4*u_p + 1/4*u_v;
b_p_P_1 = P_1(1);
b_d_P_1 = P_1(2);

trace_A_P_1_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_1 - s*omega*b_d_P_1);
trace_A_P_1 = 0;
discriminant_A_P_1_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_1 - 8*c^2*s*omega*b_d_P_1 + c^2*b_p_P_1^2 + 2*c*s*omega*b_p_P_1*b_d_P_1 + s^2*omega^2*b_d_P_1^2 - 4);
discriminant_A_P_1_2 = 4*(c^2 - s^2)^2 - 4;
discriminant_A_P_1 = 0;
lambda_P_1_1 = 0.5*(trace_A_P_1 + sqrt(discriminant_A_P_1));
lambda_P_1_2 = 0.5*(trace_A_P_1 - sqrt(discriminant_A_P_1));

% disp(simplify(trace_A_P_1_1 - trace_A_P_1))
% disp(simplify(discriminant_A_P_1_1 - discriminant_A_P_1_2))

% point 2
P_2 = P_origin + 3/2*u_p + 3/2*u_v;
b_p_P_2 = P_2(1);
b_d_P_2 = P_2(2);

trace_A_P_2_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_2 - s*omega*b_d_P_2);
trace_A_P_2 = 0;
discriminant_A_P_2_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_2 - 8*c^2*s*omega*b_d_P_2 + c^2*b_p_P_2^2 + 2*c*s*omega*b_p_P_2*b_d_P_2 + s^2*omega^2*b_d_P_2^2 - 4);
discriminant_A_P_2_2 = - 16*(c^2 - s^2)^2 - 4;
discriminant_A_P_2 = -20;
lambda_P_2_1 = 0.5*(trace_A_P_2 + sqrt(discriminant_A_P_2));
lambda_P_2_2 = 0.5*(trace_A_P_2 - sqrt(discriminant_A_P_2));

disp(simplify(trace_A_P_2_1 - trace_A_P_2))
disp(simplify(discriminant_A_P_2_1 - discriminant_A_P_2_2))

% point 3
P_3 = P_origin + 3/2*u_p - 1/4*u_v;
b_p_P_3 = P_3(1);
b_d_P_3 = P_3(2);

trace_A_P_3_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_3 - s*omega*b_d_P_3);
trace_A_P_3_2 = 7/2*(c^2 - s^2);
trace_A_P_3 = 7/2;
discriminant_A_P_3_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_3 - 8*c^2*s*omega*b_d_P_3 + c^2*b_p_P_3^2 + 2*c*s*omega*b_p_P_3*b_d_P_3 + s^2*omega^2*b_d_P_3^2 - 4);
discriminant_A_P_3_2 = 41/4*(c^4 - 2*c^2*s^2 + s^4) - 4;
discriminant_A_P_3 = 41/4 - 4;
lambda_P_3_1 = 0.5*(trace_A_P_3 + sqrt(discriminant_A_P_3));
lambda_P_3_2 = 0.5*(trace_A_P_3 - sqrt(discriminant_A_P_3));

disp(simplify(trace_A_P_3_1 - trace_A_P_3_2))
disp(simplify(discriminant_A_P_3_1 - discriminant_A_P_3_2))

% point 4
P_4 = P_origin + 1/4*u_p - 1/4*u_v;
b_p_P_4 = P_4(1);
b_d_P_4 = P_4(2);

trace_A_P_4_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_4 - s*omega*b_d_P_4);
trace_A_P_4_2 = c^2 - s^2;
trace_A_P_4 = 1;
discriminant_A_P_4_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_4 - 8*c^2*s*omega*b_d_P_4 + c^2*b_p_P_4^2 + 2*c*s*omega*b_p_P_4*b_d_P_4 + s^2*omega^2*b_d_P_4^2 - 4);
discriminant_A_P_4_2 = 9*(c^2 - s^2)^2 - 4;
discriminant_A_P_4 = 5;
lambda_P_4_1 = 0.5*(trace_A_P_4 + sqrt(discriminant_A_P_4));
lambda_P_4_2 = 0.5*(trace_A_P_4 - sqrt(discriminant_A_P_4));

disp(simplify(trace_A_P_4_1 - trace_A_P_4_2))
disp(simplify(discriminant_A_P_4_1 - discriminant_A_P_4_2))

% point 5
P_5 = P_origin - 1/4*u_p - 1/4*u_v;
b_p_P_5 = P_5(1);
b_d_P_5 = P_5(2);

trace_A_P_5_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_5 - s*omega*b_d_P_5);
trace_A_P_5 = 0;
discriminant_A_P_5_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_5 - 8*c^2*s*omega*b_d_P_5 + c^2*b_p_P_5^2 + 2*c*s*omega*b_p_P_5*b_d_P_5 + s^2*omega^2*b_d_P_5^2 - 4);
discriminant_A_P_5_2 = 12*(c^2 - s^2)^2 - 4;
discriminant_A_P_5 = 8;
lambda_P_5_1 = 0.5*(trace_A_P_5 + sqrt(discriminant_A_P_5));
lambda_P_5_2 = 0.5*(trace_A_P_5 - sqrt(discriminant_A_P_5));

disp(simplify(discriminant_A_P_5_1 - discriminant_A_P_5_2))

% point 4
P_6 = P_origin - 1/4*u_p + 1/4*u_v;
b_p_P_6 = P_6(1);
b_d_P_6 = P_6(2);

trace_A_P_6_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_6 - s*omega*b_d_P_6);
trace_A_P_6_2 = s^2 - c^2;
trace_A_P_6 = -1;
discriminant_A_P_6_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_6 - 8*c^2*s*omega*b_d_P_6 + c^2*b_p_P_6^2 + 2*c*s*omega*b_p_P_6*b_d_P_6 + s^2*omega^2*b_d_P_6^2 - 4);
discriminant_A_P_6_2 = 9*(c^2 - s^2)^2 - 4;
discriminant_A_P_6 = 5;
lambda_P_6_1 = 0.5*(trace_A_P_6 + sqrt(discriminant_A_P_6));
lambda_P_6_2 = 0.5*(trace_A_P_6 - sqrt(discriminant_A_P_6));

disp(simplify(trace_A_P_6_1 - trace_A_P_6_2))
disp(simplify(discriminant_A_P_6_1 - discriminant_A_P_6_2))

% point 3
P_7 = P_origin - 1/4*u_p + 3/2*u_v;
b_p_P_7 = P_7(1);
b_d_P_7 = P_7(2);

trace_A_P_7_1 = simplify(2*c^2 + 2*s^2 - c*b_p_P_7 - s*omega*b_d_P_7);
trace_A_P_7_2 = -7/2*(c^2 - s^2);
trace_A_P_7 = -7/2;
discriminant_A_P_7_1 = simplify(4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p_P_7 - 8*c^2*s*omega*b_d_P_7 + c^2*b_p_P_7^2 + 2*c*s*omega*b_p_P_7*b_d_P_7 + s^2*omega^2*b_d_P_7^2 - 4);
discriminant_A_P_7_2 = 41/4*(c^4 - 2*c^2*s^2 + s^4) - 4;
discriminant_A_P_7 = 41/4 - 4;
lambda_P_7_1 = 0.5*(trace_A_P_7 + sqrt(discriminant_A_P_7));
lambda_P_7_2 = 0.5*(trace_A_P_7 - sqrt(discriminant_A_P_7));

disp(simplify(trace_A_P_7_1 - trace_A_P_7_2))
disp(simplify(discriminant_A_P_7_1 - discriminant_A_P_7_2))



% ----------------------------------------------------------------------------------------------------------------------
% look at A for point 1
% ----------------------------------------------------------------------------------------------------------------------

% look at A, because this is the optimal point
A_1 = [(- c*b_p_P_1 + c^2 + s^2),  - (c*b_d_P_1 - 2/omega*c*s); (2*c - b_p_P_1)*omega*s, (- b_d_P_1*omega*s + c^2 + s^2)];  % definition of A
A_2 = [(- c*(s^2/(2*c) + (3*c)/2) + c^2 + s^2),  - (c*(c^2/(2*omega*s) + (3*s)/(2*omega)) - 2/omega*c*s); (2*c - (s^2/(2*c) + (3*c)/2))*omega*s, (- (c^2/(2*omega*s) + (3*s)/(2*omega))*omega*s + c^2 + s^2)];  % plug in b_p and b_d
A_3 = [-1/2*(c^2 - s^2), -(c*(c^2 - s^2))/(2*omega*s); (omega*s*(c^2 - s^2))/(2*c), 1/2*(c^2 - s^2)];     % simplify
A_4 = [-1/2*(1        ), -(c*(1        ))/(2*omega*s); (omega*s*(1        ))/(2*c), 1/2*(1        )];     % apply hyperbolic identity
A_5 = [-1/2, -c/(2*omega*s); omega*s/(2*c), 1/2];     % simplify

% disp(simplify(A_1 - A_2))
% disp(simplify(A_1 - A_3))
% disp(simplify(A_4 - A_5))












