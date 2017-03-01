#! /bin/sh

echo "\nRunning the command: mpirun -np $num_slots $TA_BIN/align_tiles_split $mat_file\n"
time mpirun -np $num_slots $TA_BIN/align_tiles_split $mat_file $params_file


