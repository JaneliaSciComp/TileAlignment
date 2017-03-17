% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function merge_mat_files(target_mat_file, num_workers)
    A = [];
    b = [];
    for i=1:num_workers
        input_mat_file = [ target_mat_file(1:(length(target_mat_file)-4)) '_' num2str(i) '.mat'];
        mat = load(input_mat_file);
        A = [A   mat.A];
        b = [b ; mat.b];
        disp(['    i=' num2str(i) ' size(A)=' num2str(size(A)) ' size(b)=' num2str(size(b))]);
    end
    save(target_mat_file, 'A', 'b', '-v7.3')
    disp(['produced file ' target_mat_file ])
