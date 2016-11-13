#! /bin/sh

echo "\nRunning the command: mpirun -np $num_slots $TA_HOME/align_multi $TA_DATA/$mat_file\n"
time mpirun -np $num_slots $TA_HOME/align_multi $TA_DATA/$mat_file $params_file


