function count_tiles(input_mat_file)

    mat = load(input_mat_file);
    num_col = size(mat.A, 2);
    num_tiles = num_col/6;
    disp(['num_tiles=' num2str(num_tiles)]);
