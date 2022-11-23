function run()
	clc; tic;

    % transition to first energy level with different values of A and w
	A_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.5];
	w_list = [0.9, 0.95, 0.98, 0.99, 1.01, 1.02, 1.05, 1.1];
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
    fprintf("Part 1 done\n");

    % two-photon transition
	main(0.01, 0.49, 40*pi);
	fprintf("A = %.3f, w = %.3f done\n", 0.01, 0.49);

    fprintf("Part 2 done\n");

    toc;

    % omega ~ 2 for 0->2 and 0->4
    w_list = [1.7, 1.8, 1.9, 1.95, 1.98, 2.02, 2.05, 2.1, 2.2, 2.3];
    state_num_list = [2,4];
    t_final = 120 * pi;

    w_index = 1 : length(w_list);
    state_num_index = 1 : length(state_num_list);

    parfor i1 = w_index
        for i2 = state_num_index
	        main(0.01, w_list(i1), t_final, state_num_list(i2));
            fprintf("A = 0.01, w = %.3f, State 0-%i done\n", w_list(i1), state_num_list(i2));
        end
    end

    fprintf("All done\n");
	toc;


    
end

