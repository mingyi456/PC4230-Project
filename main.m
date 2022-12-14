function [time_axis, Prob, Prob_theo] = main(A, w, t_final, target_state_num)
	hold on;
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

    fig_prob_t = figure(1);
    fig_prob_t.Position = [10 10 900 900];
	time_axis = T * dt / pi;
	plot(time_axis, Prob);

	v_matrix = calc_v_matrix(A, 20);
	w_nm = target_state_num;
	
	if target_state_num == 1
		v_nm = v_matrix(1, target_state_num + 1);
		Prob_theo_prop = (abs(v_nm) ^ 2) / 4;
		major_term = (exp(1i * (w_nm - w) * T * dt) - 1) ./ (w_nm - w);
		minor_term = (exp(1i * (w_nm + w) * T * dt) - 1) ./ (w_nm + w);
		Prob_theo = Prob_theo_prop * abs(major_term + minor_term) .^ 2;
		plot(time_axis, Prob_theo);
		legend(["Simulated", "Theoretical"], FontSize = 15);

    elseif target_state_num == 4
		v_sum = v_matrix(target_state_num + 1, :) * (v_matrix(:, 1) ./ (v_matrix(:, 1) - w));
		Prob_theo_prop = (abs(v_sum) ^ 2);
		sin_arg = (w_nm - 2 * w) * (T * dt) / 2;
		Prob_theo = Prob_theo_prop * (sin(sin_arg) .^ 2) / ((w_nm - 2 * w) ^ 2);
		Prob_theo = Prob_theo * (max(Prob) / max(Prob_theo));
		plot(time_axis, Prob_theo);
        legend(["Simulated", "Theoretical"], FontSize = 15);

    else
		Prob_theo = zeros(1);
		legend(["Simulated"], FontSize=15);
	end

	title(sprintf("A = %.3f, w = %.3f for state 0\\rightarrow%i", A, w, target_state_num), FontSize=18);
	ax = gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    xlabel('t/\pi', FontSize=18);
	xlim([0,T(end)*dt/pi])
    ylabel(sprintf('P_{0\\rightarrow%i}', target_state_num), FontSize=18);
    % ylim([0,1.1*max(Prob_theo)])
    legend('boxoff')
    file_name = sprintf("A=%.3f, w=%.3f, t_final=%gpi, State=0-%i", A, w, t_final/pi, target_state_num);
	saveas(gcf, sprintf(".\\final_plots\\%s.png", file_name));
	close;
	
% 	plot(T * dt, freq(Prob));
% 	title(sprintf("A = %.3f, w = %.3f", A, w));
% 	[peak_heights , peak_positions] = findpeaks(f, T * dt)

end

function f = freq(Prob)
	f = fftshift(abs(fft(Prob)));
end

