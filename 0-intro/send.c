#include <stdio.h>
#include "mpi.h"
#include <string.h>

int main(int argc, char * argv[]){

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // printf("size: %d rank: %d\n", size, rank);

    char buf[100];
    sprintf(buf, "Hello from %d", rank);
    MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, (rank + 1) % size, 0, MPI_COMM_WORLD);
    MPI_Recv(buf, 100, MPI_CHAR, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("#%d: '%s'\n", rank, buf);

    MPI_Finalize();

    return 0;
}