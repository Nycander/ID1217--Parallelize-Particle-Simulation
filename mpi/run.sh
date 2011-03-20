#!/bin/sh

N_PROC=$((2*$1+1))
mpiexec -np $N_PROC ./mpi $@