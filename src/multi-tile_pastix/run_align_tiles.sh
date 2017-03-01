#! /bin/sh

echo "\nRunning the command: mpirun -np $num_slots $TA_BIN/align_tiles $mat_file\n"
time mpirun -np $num_slots  -machinefile /tmp/*/machines $TA_HOME/align_tiles $mat_file $params_file


