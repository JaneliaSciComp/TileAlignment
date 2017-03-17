#! /bin/bash

export I_MPI_FABRICS=shm:tcp
export I_MPI_FALLBACK=disable
qsub -v num_slots=$1,mat_file=$2,params_file=$3 -A flyTEM -R yes -N ats$1 -pe impi2 $1 -now n -V $TA_BIN/run_align_multi.sh



