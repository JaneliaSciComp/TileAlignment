% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function split_mat_file(input_mat_file, num_workers)

    mat = load(input_mat_file);
    num_col = size(mat.A, 2);
    bmax = max(abs(mat.b));
    tiles_per_worker = round(num_col/6./num_workers);
    disp(['tiles_per_worker=' num2str(tiles_per_worker)]);
    for i=1:num_workers
        col_min = 1 + 6*(i-1)*tiles_per_worker;
        if i < num_workers
            col_max = col_min   + 6*tiles_per_worker-1;
        else
            col_max = size(mat.A, 2);
        end
        disp(['    i=' num2str(i) ' col_min=' num2str(col_min) ' col_max=' num2str(col_max)]);
        prefixName = input_mat_file(1:(length(input_mat_file)-4));
        output_mat_file = strcat(prefixName,'_', num2str(i), '.mat');
        A = mat.A(:, col_min:col_max);
        b = mat.b(col_min:col_max);
        disp(['    i=' num2str(i) ' size(A)=' num2str(size(A)) ' size(b)=' num2str(size(b))]);
        save(output_mat_file, 'A', 'b', 'bmax', '-v7.3')
    end
