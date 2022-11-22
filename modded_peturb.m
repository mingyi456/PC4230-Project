function Prob = modded_peturb(A, w, t_final, target_state_num, peturb_type)
	clc; hold on;
	x_min = -80; x_max = 80;
	x_range = x_max - x_min;
	num_points = 2048;
	X = x_min + x_range * (0 : num_points - 1) / num_points; 
	P = (2 * pi / x_range) * [0 : num_points / 2 - 1,-num_points / 2 : -1];
	
	gaussian_state_centre = 0;
	gaussian_state_width = 1;
	gaussian_state = exp(-(X - gaussian_state_centre) .^2 / (2 * gaussian_state_width ^ 2));

	if ~exist('t_final','var')
      	t_final = 400;
	end
	if ~exist('target_state_num','var')
      	target_state_num = 1;
	end

    if ~exist('peturb_type','var')
        peturb_type = 'none';
    end

	initial_state = normalize(hermiteH(0 , X) .* gaussian_state, "norm");
	target_state = normalize(hermiteH(target_state_num , X) .* gaussian_state, "norm");

	dt = 0.005;
	num_steps = ceil(t_final / dt); T = 1 : num_steps;
	
	Prob = zeros(1, num_steps);


	U_t = exp(-1i * dt * (P .^ 2 / 2));
	curr_state = initial_state;
	
    if strcmp(peturb_type,'none')
        peturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(peturb_type,'sawtooth')
        peturb_term = @(x) A * sawtooth(X-pi) * sawtooth(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(peturb_type,'triangle')
        peturb_term = @(x) A * sawtooth(X+pi/2,1/2) * -sawtooth(w * x * dt,1/2) + (X .^ 2 / 2);
    elseif strcmp(peturb_type,'square')
        peturb_term = @(x) A  * square(X) * square(w * x * dt) + (X .^ 2 / 2);
    elseif strcmp(peturb_type,'quadratic')
        peturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2) + (X .^ 4 / 8);
    else
        peturb_term = @(x) A * sin(X) * cos(w * x * dt) + (X .^ 2 / 2);
    end

	for i1 = T
    	U_v = exp(((-1i/2) * dt) * peturb_term(i1));

		curr_state = U_v .* curr_state;
		curr_state = U_t .* fft(curr_state);
		curr_state = U_v .* ifft(curr_state);
	
		dot_pdt = target_state * transpose(curr_state);
		Prob(i1) = dot_pdt * dot_pdt';
	
	end

    plot(T * dt/pi, Prob);
    title(sprintf("A = %.3f, w = %.3f for state 0\\rightarrow%i, %s", A, w, target_state_num, peturb_type), FontSize=18);
    xlabel('t/\pi', FontSize=18);
    xlim([0,T(end)*dt/pi])
    ylabel(sprintf('P_{0\\rightarrow%i}', target_state_num), FontSize=18);
    file_name = sprintf("A=%.3f, w=%.3f, t_final=%g, State=0-%i, %s", A, w, t_final, target_state_num, peturb_type);
    saveas(gcf, sprintf(".\\plots\\%s.svg", file_name));
	% close;
	
% 	plot(T * dt, freq(Prob));
% 	title(sprintf("A = %.3f, w = %.3f", A, w));
% 	[peak_heights , peak_positions] = findpeaks(f, T * dt)

end

function f = freq(Prob)
	f = fftshift(abs(fft(Prob)));
end

