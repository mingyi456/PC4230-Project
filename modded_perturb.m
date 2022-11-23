function [time_axis, Prob] = modded_perturb(A, w, t_final, target_state_num, perturb_type)
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
    if ~exist('perturb_type','var')
        perturb_type = 'none';
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

    if strcmp(perturb_type,'none')
        perturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(perturb_type,'sawtooth')
        perturb_term = @(x) A * sawtooth(X-pi) * sawtooth(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(perturb_type,'triangle')
        perturb_term = @(x) A * sawtooth(X+pi/2,1/2) * -sawtooth(w * x * dt,1/2) + (X .^ 2 / 2);
    elseif strcmp(perturb_type,'square')
        perturb_term = @(x) A  * square(X) * square(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(perturb_type,'quartic')
        perturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2) + (X .^ 4 / 8);
    else
        perturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2);
    end

	for i1 = T
		U_v = exp((-1i/2) * dt * perturb_term(i1));
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

% 	v_matrix = calc_v_matrix(A,10);
% 	w_nm = target_state_num;
	
	title(sprintf("A = %.3f, w = %.3f for state 0\\rightarrow%i (%s)", A, w, target_state_num), perturb_type, FontSize=18);
	ax = gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    xlabel('t/\pi', FontSize=18);
	xlim([0,T(end)*dt/pi])
    ylabel(sprintf('P_{0\\rightarrow%i}', target_state_num), FontSize=18);
    legend(["Simulated"],'boxoff')
    file_name = sprintf("A=%.3f, w=%.3f, t_final=%gpi, State=0-%i, %s", A, w, t_final/pi, target_state_num, perturb_type);
	saveas(gcf, sprintf(".\\final_plots\\%s.png", file_name));
	close;
	
% 	plot(T * dt, freq(Prob));
% 	title(sprintf("A = %.3f, w = %.3f", A, w));
% 	[peak_heights , peak_positions] = findpeaks(f, T * dt)

end

function f = freq(Prob)
	f = fftshift(abs(fft(Prob)));
end

