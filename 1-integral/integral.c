#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double f(double x) { return sqrt(4 - x * x); }

int main(int argc, char *argv[]) {

    int rank, size;
    unsigned int N = atoi(argv[1]);
    double h = 2. / N;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start_time, end_time;
    if (size == 1) {
        start_time = MPI_Wtime();

        double res = 0.;
        for (unsigned int i = 0; i < N; i++) {
            res += (f(h * (i + 1)) + f(h * i)) * h / 2;
        }
        // printf("%.6f", res);
        end_time = MPI_Wtime();
        printf("%f\n", end_time - start_time);
    } else if (rank == size - 1) {
        start_time = MPI_Wtime();
        double res = 0.;
        double res_buff[size - 1];
        unsigned int interval_size = ceil((double) N / size);
        unsigned int last_interval_size = N - (size - 1) * interval_size;
        for (unsigned int i = N - 1; i > N - last_interval_size - 1; i--) {
            res += (f(h * (i + 1)) + f(h * i)) * h / 2;
        }
        for (int i = 0; i < size - 1; i++) {
            MPI_Recv(&res_buff[i], 1, MPI_DOUBLE, i, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            res += res_buff[i];
        }
        // printf("%.6f\n", res);
        end_time = MPI_Wtime();
        printf("%f\n", end_time - start_time);
    } else {
        double res = 0.;
        unsigned int interval_size = ceil((double) N / size);
        for (unsigned int i = rank * interval_size;
             i < (rank + 1) * interval_size; i++) {
            res += (f(h * (i + 1)) + f(h * i)) * h / 2;
        }
        MPI_Send(&res, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}