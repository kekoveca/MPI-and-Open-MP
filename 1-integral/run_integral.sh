#!/bin/zsh
# Loop from 1 to 28
for i in {1..28}; do
    echo "\nRunning mpirun with -np $i"
    gtime -f "%e" -o times.txt -a mpirun --map-by :OVERSUBSCRIBE -np $i ./mpi_integral 1000000000000
done

echo "All runs completed. Check times.txt for results."