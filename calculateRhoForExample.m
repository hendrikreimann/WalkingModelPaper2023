cadence_to_use = 120;
minute_to_second = 1/60;
bpm = cadence_to_use;
metronome_frequency = bpm * 0.5 * minute_to_second; % two beats per metronome cycle
stride_time = metronome_frequency.^(-1);
step_time = stride_time * 0.5;

% set parameters
z_c = 1;
g = 9.81;

% calculate dependent variables
T_total = 1/metronome_frequency * 0.5;
T_step = 1/(2*metronome_frequency);
omega = sqrt(g/z_c);
c = cosh(omega * T_step * 0.5);
s = sinh(omega * T_step * 0.5);

b_p = c + s;
b_d = 1/omega * (c + s);

A_here = ...
    [ ...
      (- c*b_p + c^2 + s^2),  - (c*b_d - 2/omega*c*s); ...
      (2*c - b_p)*omega*s, (- b_d*omega*s + c^2 + s^2) ...
    ];

% Eigendecomposition to get spectral norm rho(M) (https://en.wikipedia.org/wiki/Convergent_matrix)
D = eig(A_here);
rho = max(abs(D));





















