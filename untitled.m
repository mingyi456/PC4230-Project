hold on;
x_min = -80; x_max = 80;
x_range = x_max - x_min;
num_points = 2048;
X = x_min + x_range * (0 : num_points - 1) / num_points; 
P = (2 * pi / x_range) * [0 : num_points / 2 - 1,-num_points / 2 : -1];

t_final = 18;

dt = 0.005;
num_steps = ceil(t_final / dt); T = 1 : num_steps;

Prob = zeros(1, num_steps);
% U_v_0 = exp((-1i * dt / 2 ) * (X .^ 2 / 2));
U_t = exp(-1i * dt * (P .^ 2 / 2));

gaussian_state_centre = 0;
gaussian_state_width = 1;
gaussian_state = exp(-(X - gaussian_state_centre) .^2 / (2 * gaussian_state_width ^ 2));

hermite_polynorm = hermiteH(0 , X);
initial_state = normalize(hermite_polynorm .* gaussian_state, "norm");

excited_state_1 = normalize(hermiteH(1, X) .* gaussian_state, "norm");
excited_state_2 = normalize(hermiteH(2, X) .* gaussian_state, "norm");

assert(conj(initial_state) * initial_state' == 1, "Initial state not normalised");


A = 0.01;
w = 0.48;


curr_state = initial_state;

for i1 = T
	U_v = exp(((-1i/2) * dt) * ((A * sin(X) * cos(w * i1 * dt)) + (X .^ 2 / 2)));

	curr_state = U_v .* curr_state;
	curr_state = U_t .* fft(curr_state);
	curr_state = U_v .* ifft(curr_state);

	dot_pdt = conj(excited_state_1) * curr_state';
	Prob(i1) = conj(dot_pdt) * dot_pdt;

end

plot(T * dt, Prob);
title(sprintf("A = %.3f, w = %.3f", A, w));
file_name = sprintf("A=%.3f, w=%.3f, t_final=%g", A, w, t_final);
saveas(gcf, sprintf(".\\plots\\%s.svg", file_name));