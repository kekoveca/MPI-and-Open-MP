#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);               // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the current process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processes

    if (size < 2) {
        if (rank == 0) {
            printf("This program requires at least two processes.\n");
        }
        MPI_Finalize();
        return 0;
    }

    const int tag = 0;
    int max_bytes = 10000000;
    double start_time, end_time;
    MPI_Status status;

    if (rank == 0) {
        printf("Buffer,Time\n");
    }

    for (int buffer_size = 1; buffer_size <= max_bytes; buffer_size *= 2) {

        if (rank == 0) {
            char* buffer = (char*) malloc(buffer_size * sizeof(char));
            start_time = MPI_Wtime();
            MPI_Send(buffer, buffer_size, MPI_BYTE, 1, tag, MPI_COMM_WORLD);
            MPI_Recv(buffer, buffer_size, MPI_BYTE, 1, tag, MPI_COMM_WORLD, &status);
            end_time = MPI_Wtime();

            double total_time = (end_time - start_time) / 2;

            printf("%d,%.6f\n", buffer_size, total_time);
            free(buffer);

        } else if (rank == 1) {
            char* buffer = (char*) malloc(buffer_size * sizeof(char));
            MPI_Recv(buffer, buffer_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Send(buffer, buffer_size, MPI_BYTE, 0, tag, MPI_COMM_WORLD);
            free(buffer);
        }
    }

    MPI_Finalize(); 
    return 0;
}