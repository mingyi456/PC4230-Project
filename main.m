function [time_axis, Prob, Prob_theo] = main(A, w, t_final, target_state_num) 
	clc; hold on;
	x_min = -80; x_max = 80;
	x_range = x_max - x_min;
	num_points = 2048;
	X = x_min + x_range * (0 : num_points - 1) / num_points; 
	P = (2 * pi / x_range) * [0 : num_points / 2 - 1,-num_points / 2 : -1];
	
	gaussian_state_centre = 0;
	gaussian_state_width = 1;
	gaussian_state = exp(-(X - gaussian_state_centre) .^2 / (2 * gaussian_state_width ^ 2));

	if ~exist('target_state_num','var')
      	target_state_num = 1;
	end

	initial_state = normalize(hermiteH(0 , X) .* gaussian_state, "norm");
	target_state = normalize(hermiteH(target_state_num , X) .* gaussian_state, "norm");

	if ~exist('t_final','var')
      	t_final = 400;
	end

	dt = 0.005;
	num_steps = ceil(t_final / dt); T = 1 : num_steps;
	
	U_t = exp(-1i * dt * (P .^ 2 / 2));
	curr_state = initial_state;
	Prob = zeros(1, num_steps);

	for i1 = T
		U_v = exp(((-1i/2) * dt) * ((A * sin(X) * cos(w * i1 * dt)) + (X .^ 2 / 2)));
		curr_state = U_v .* curr_state;
		curr_state = U_t .* fft(curr_state);
		curr_state = U_v .* ifft(curr_state);
		dot_pdt = target_state * transpose(curr_state);
		Prob(i1) = dot_pdt * dot_pdt';
	
	end

	v_matrix = calc_v_matrix(A);
	w_nm = target_state_num;
	
	if rem(target_state_num, 2) ~= 0
		v_nm = v_matrix(1, target_state_num + 1);
		Prob_theo_prop = 2 * (abs(v_nm) ^ 2);
		Prob_theo_major = Prob_theo_prop * (sin((w - w_nm) * T * dt / 2).^2) / ((w - w_nm) ^ 2);
		Prob_theo_minor = Prob_theo_prop * (sin((w + w_nm) * T * dt / 2).^2) / ((w + w_nm) ^ 2);
		Prob_theo = Prob_theo_major + Prob_theo_minor;
	else
		v_sum = v_matrix(target_state_num + 1, :) * (v_matrix(:, 1) ./ (v_matrix(:, 1) - w));
		Prob_theo_prop = (abs(v_sum * 40) ^ 2);
		Prob_theo = Prob_theo_prop * (sin((w_nm - 2 * w) * T * dt / 2) .^ 2) / ((w_nm - 2 * w) ^ 2);
	end

	time_axis = T * dt / pi;
			
	plot(time_axis, Prob);
	plot(time_axis, Prob_theo);
	legend(["Simulated Prob", "Theoretical Prob"]);
	title(sprintf("A = %.3f, w = %.3f for state 0\\rightarrow%i", A, w, target_state_num), FontSize=18);
	xlabel('t/\pi', FontSize=18);
	xlim([0,T(end)*dt/pi])
    ylabel(sprintf('P_{0\\rightarrow%i}', target_state_num), FontSize=18);
	file_name = sprintf("A=%.3f, w=%.3f, t_final=%gpi, State=0-%i", A, w, t_final/pi, target_state_num);
	saveas(gcf, sprintf(".\\plots\\%s.svg", file_name));
	% close;
	
% 	plot(T * dt, freq(Prob));
% 	title(sprintf("A = %.3f, w = %.3f", A, w));
% 	[peak_heights , peak_positions] = findpeaks(f, T * dt)

end

function f = freq(Prob)
	f = fftshift(abs(fft(Prob)));
end

