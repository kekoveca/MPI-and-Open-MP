#!/bin/bash
#PBS -l walltime=00:10:00,nodes=7:ppn=4
#PBS -N mpi_life_string
#PBS -q batch

cd $PBS_O_WORKDIR
for i in {1..28}; do
    mpirun --hostfile $PBS_NODEFILE -np $i ./mpi_life_string p46gun_big.cfg
done