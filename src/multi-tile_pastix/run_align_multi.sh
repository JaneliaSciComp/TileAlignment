#! /bin/sh

echo "\nRunning the command: mpirun -np $num_slots $TA_BIN/align_multi $mat_file\n"
time mpirun -np $num_slots  -machinefile /tmp/*/machines $TA_BIN/align_multi $mat_file $params_file


