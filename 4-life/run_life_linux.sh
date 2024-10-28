#!/bin/bash
# Loop from 1 to 12
for i in {1..12}; do
    echo "\nRunning mpirun with -np $i"
    \time -f "%e" -o times.txt -a mpirun --map-by :OVERSUBSCRIBE -np $i ./life_mpi p46gun_big.cfg
done

echo "All runs completed. Check times.txt for results."