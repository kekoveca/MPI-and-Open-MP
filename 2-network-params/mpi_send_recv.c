#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        if (rank == 0) {
            printf("This program requires at least two processes.\n");
        }
        MPI_Finalize();
        return 0;
    }

    int max_size = 1e6;
    int num_repeats = 1e5;

    for (double msg_size = 1; msg_size <= max_size; msg_size *= 10) {

        char *msg_buffer = (char *) calloc(msg_size, sizeof(char));
        double time = MPI_Wtime();
        for (int i = 0; i < num_repeats; ++i) {
            if (rank == 0) {
                MPI_Send(msg_buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(msg_buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (rank == 1) {
                MPI_Recv(msg_buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(msg_buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            time = MPI_Wtime() - time;
            printf("%f,%f\n", msg_size, 0.5 * time / num_repeats * 1e6);
        }

        free(msg_buffer);
    }

    MPI_Finalize();
    return 0;
}