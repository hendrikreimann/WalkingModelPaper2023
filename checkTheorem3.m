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
M = [(c^2 + s^2 - c*b_p), (2*s/omega - b_d)*c; (2*c - b_p)*s*omega, (c^2 + s^2 - s*omega*b_d)];

% ----------------------------------------------------------------------------------------------------------------------
% determinant
% ----------------------------------------------------------------------------------------------------------------------
det_M_1 = M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1); % definition of the determinant
det_M_2 = (c^2 + s^2 - c*b_p)*(c^2 + s^2 - s*omega*b_d) - (2*s/omega - b_d)*c*(2*c - b_p)*s*omega; % plug in definition from above
det_M_3 = c^4 + 2*c^2*s^2 - b_p*c^3 + s^4 - c*s^2*b_p - c^2*s*omega*b_d - s^3*omega*b_d + c*s*omega*b_p*b_d - 4*c^2*s^2 + 2*c^2*s*omega*b_d + 2*c*s^2*b_p - c*s*omega*b_d*b_p; % expand
det_M_4 = (c^2 - s^2)^2 - (c^2 - s^2)*c*b_p + (c^2 - s^2)*s*omega*b_d; % sort
det_M_5 = 1             - 1          *c*b_p + 1*          s*omega*b_d; % apply hyperbolic identity

% verify transformations
% disp(simplify(det_M_1 - det_M_2))
% disp(simplify(det_M_1 - det_M_3))
% disp(simplify(det_M_1 - det_M_4))

% ----------------------------------------------------------------------------------------------------------------------
% trace
% ----------------------------------------------------------------------------------------------------------------------

trace_M_1 = M(1, 1) + M(2, 2);
trace_M_2 = (c^2 + s^2 - c*b_p) + (c^2 + s^2 - s*omega*b_d);    % plug in definition from above
trace_M_3 = c^2 + s^2 - c*b_p + c^2 + s^2 - s*omega*b_d;        % expand
trace_M_4 = 2*(c^2 + s^2) - (c*b_p + s*omega*b_d);              % sort
trace_M   = 2*c^2 + 2*s^2 - c*b_p - s*omega*b_d;

% verify transformations
% disp(simplify(trace_M_1 - trace_M_2))
% disp(simplify(trace_M_1 - trace_M_3))
% disp(simplify(trace_M_1 - trace_M_4))
% disp(simplify(trace_M_1 - trace_M))

% ----------------------------------------------------------------------------------------------------------------------
% Discriminant
% ----------------------------------------------------------------------------------------------------------------------

discriminant_M_1 = trace_M_4^2 - 4*det_M_5;                                                     % definition of the discriminant
discriminant_M_2 = (2*(c^2 + s^2) - (c*b_p + s*omega*b_d))^2 - 4*(1 - c*b_p + s*omega*b_d);     % plug in equations from above
discriminant_M_3 = (2*c^2 + 2*s^2 - c*b_p - s*omega*b_d)*(2*c^2 + 2*s^2 - c*b_p - s*omega*b_d) - 4 + 4*c*b_p - 4*s*omega*b_d;     % expand
discriminant_M_4 = 4*c^4 + 4*s^2*c^2 - 2*c^3*b_p - 2*c^2*s*omega*b_d ...
                   + 4*c^2*s^2 + 4*s^4 - 2*c*s^2*b_p - 2*s^3*omega*b_d ...
                   - 2*c^3*b_p - 2*c*s^2*b_p + c^2*b_p^2 + c*s*omega*b_p*b_d ...
                   - 2*c^2*s*omega*b_d - 2*s^3*omega*b_d + c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 ...
                   - 4 + 4*c*b_p - 4*s*omega*b_d;     % expand
discriminant_M_5 = 4*(c^2 + s^2)^2 - 4*(c^2 + s^2 - 1)*c*b_p - 4*(c^2 + s^2 + 1)*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;     % sort
discriminant_M_6 = 4*(c^2 + s^2)^2 - 4*(2*s^2 + 1 - 1)*c*b_p - 4*(2*c^2 - 1 + 1)*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;     % apply hyperbolic identities
discriminant_M_7   = 4*(c^2 + s^2)^2 - 8*s^2*c*b_p - 8*c^2*s*omega*b_d + (c*b_p + s*omega*b_d)^2 - 4;     % simplify
discriminant_M     = 4*c^4 + 4*2*c^2*s^2 + 4*s^4 - 8*s^2*c*b_p - 8*c^2*s*omega*b_d + c^2*b_p^2 + 2*c*s*omega*b_p*b_d + s^2*omega^2*b_d^2 - 4;

% verify transformations
% disp(simplify(discriminant_M_1 - discriminant_M_2))
% disp(simplify(discriminant_M_1 - discriminant_M_3))
% disp(simplify(discriminant_M_1 - discriminant_M_4))
% disp(simplify(discriminant_M_1 - discriminant_M_5))
% disp(simplify(discriminant_M_6 - discriminant_M_7))
% disp(simplify(discriminant_M_6 - discriminant_M))

% ----------------------------------------------------------------------------------------------------------------------
% Eigenvalue - discriminant positive
% ----------------------------------------------------------------------------------------------------------------------

% define Eigenvalue equation
eq_eigenvalues_dpos_1 = trace_M^4 - 2*trace_M^2*discriminant_M - 8*trace_M^2  + (discriminant_M-4)^2;
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
eq_eigenvalues_dneg_1 = 1/4*trace_M^2 + 1/4*(-discriminant_M) - 1;
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
disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_2))
disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_3))
disp(simplify(eq_eigenvalues_dneg_1 - eq_eigenvalues_dneg_4))


















