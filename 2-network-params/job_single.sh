#!/bin/bash
#PBS -l walltime=00:10:00,nodes=1:ppn=2
#PBS -N mpi_single
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 2 ./mpi_send_recv