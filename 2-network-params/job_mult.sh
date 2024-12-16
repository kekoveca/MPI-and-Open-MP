#!/bin/bash
#PBS -l walltime=00:10:00,nodes=2:ppn=1
#PBS -N mpi_mult
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 2 ./mpi_send_recv