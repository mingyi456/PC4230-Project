function run()
	clc; tic;

    % transition to first energy level with different values of A and w
	A_list = [0.01, 0.05, 0.1, 0.2, 0.5];
	w_list = [0.9,1.1];
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
    fprint("Part 1 done");

    % two-photon transition
	main(0.01, 0.49, 40*pi);
	fprintf("A = %.3f, w = %.3f done\n", 0.01, 0.49);

    fprint("Part 2 done");
    toc;

    % write v_matrix to csv
    v_matrix = calc_v_matrix(0.01, 100);
    writematrix(v_matrix,'v_matrix, A=0.01.csv');

    fprint("Part 3 done");
    toc;

    % omega ~ 2 for 0->2 and 0->4
    w_list = [1.9, 1.7];
    state_num_list = [2,4];
    t_final = 120 * pi;

    w_index = 1 : length(w_list);
    state_num_index = 1 : length(state_num_list);

    parfor i1 = w_index
        for i2 = state_num_index
	        main(A, w, t_final, state_num_list(i2));
            fprintf("A = 0.01, w = %.3f, State 0-%i done\n", w_list(i1), state_num_list(i2));
        end
    end

    fprintf("All done\n");
	toc;


    
end

