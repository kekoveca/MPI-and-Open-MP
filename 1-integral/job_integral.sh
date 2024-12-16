#!/bin/bash
#PBS -l walltime=00:01:00,nodes=7:ppn=4
#PBS -N mpi_integral
#PBS -q batch

cd $PBS_O_WORKDIR
for i in {1..28}; do
    mpirun --hostfile $PBS_NODEFILE -np $i ./mpi_integral 1000000000000
done