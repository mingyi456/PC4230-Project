function run()
	clc; tic;
	A_list = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.12];

	w_list = [0.45, 0.48, 0.49, 0.51, 0.52];

	t_final = 120 * pi;

	a_index = 1 : length(A_list);
	w_index = 1 : length(w_list);

	parfor i1 = w_index
		for i2 = a_index
			main(A_list(i2), w_list(i1), t_final);
			fprintf("A = %.3f, w = %.3f done\n", A_list(i2), w_list(i1));
		end
	
	end

	toc;
	fprintf("All done\n");

end

