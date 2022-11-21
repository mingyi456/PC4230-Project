clear; clc; tic;
x_min = -80; x_max = 80;
x_range = x_max - x_min;
num_points = 2048;
X = x_min + x_range * (0 : num_points - 1) / num_points; 

gaussian_state_centre = 0;
gaussian_state_width = 1;
gaussian_state = exp(-(X - gaussian_state_centre) .^2 / (2 * gaussian_state_width ^ 2));

size = 40;
A = 1;
V = A * sin(X);

states = zeros(size, num_points);
parfor i1 = 1 : size
	states(i1, :) = normalize(hermiteH(i1 -1, X) .* gaussian_state, "norm");
end

v_matrix = zeros(size, size);

parfor m = 1 : size
	row = zeros(1, size);
	state_m = states(m, :);
	for n = m : size
		if rem(m - n, 2) ~= 0 
			state_n = states(n, :);
			V_mn = state_m * transpose(V .* state_n);
			row(n) = V_mn;
		end
	end
	v_matrix(m, :) = row;
end

v_matrix = v_matrix + transpose(triu(v_matrix, 1));

toc;
v_matrix
