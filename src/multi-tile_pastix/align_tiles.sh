#!/bin/bash

export I_MPI_FABRICS=shm:ofa
export I_MPI_FALLBACK=disable
qsub -v num_slots=$1,mat_file=$2,params_file=$3 -A flyTEM -R yes -N at -pe impi $1 -now n -V $TA_BIN/run_align_tiles.sh



