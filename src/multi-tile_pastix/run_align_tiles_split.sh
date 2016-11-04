#! /bin/sh

echo "\nRunning the command: mpirun -np $num_slots $PASTIX_HOME/align_tiles_split_pastix $PASTIX_DATA/$mat_file\n"
time mpirun -np $num_slots $PASTIX_HOME/align_tiles_split_pastix $PASTIX_DATA/$mat_file $params_file


