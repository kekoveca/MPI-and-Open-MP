/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, 2016
 */

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define ind(i, j) (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx))

typedef struct {
    int nx, ny;
    int *u0;
    int *u1;
    int steps;
    int save_steps;

    /* MPI */
    int rank;
    int size;
    int start[2];
    int stop[2];
    int coords[2];
    int dims[2];

    MPI_Comm comm_cart;
    MPI_Comm comm_dims[2];

    MPI_Datatype block_t;
    MPI_Datatype row_t;

    MPI_Datatype column;
    MPI_Datatype row;

    double start_time, end_time;
} life_t;

void life_init(const char *path, life_t *l);
void life_free(life_t *l);
void life_step(life_t *l);
void life_save_vtk(const char *path, life_t *l);
void decomposition(const int n, const int p, const int k, int *start, int *stop);

void life_collect(life_t *l);
void life_collect_1d(life_t *l, MPI_Comm comm, MPI_Datatype block_t, int ind_send, int ind_recv);
void exchange_rows(life_t *l);
void exchange_columns(life_t *l);
void exchange_corners(life_t *l);
void life_exchange(life_t *l);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    if (argc != 2) {
        printf("Usage: %s input file.\n", argv[0]);
        return 0;
    }
    life_t l;
    life_init(argv[1], &l);

    int i;
    char buf[100];
    if (l.rank == l.size - 1)
        l.start_time = MPI_Wtime();
    for (i = 0; i < l.steps; i++) {
        if (i % l.save_steps == 0) {
            life_collect(&l);
            if (l.rank == l.size - 1) {
                sprintf(buf, "vtk/life_%06d.vtk", i);
                // printf("Saving step %d to '%s'.\n", i, buf);
                life_save_vtk(buf, &l);
            }
        }
        life_exchange(&l);
        life_step(&l);
    }

    if (l.rank == l.size - 1) {
        l.end_time = MPI_Wtime();
        printf("%f\n", l.end_time - l.start_time);
    }

    life_free(&l);
    MPI_Finalize();
    return 0;
}

void life_collect(life_t *l) {
    life_collect_1d(l, l->comm_dims[0], l->block_t, ind(l->start[0], l->start[1]), ind(0, l->start[1]));
    life_collect_1d(l, l->comm_dims[1], l->row_t, ind(0, l->start[1]), ind(0, 0));
}

void life_collect_1d(life_t *l, MPI_Comm comm, MPI_Datatype block_t, int ind_send, int ind_recv) {
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Gather(l->u0 + ind_send, 1, block_t, l->u0 + ind_recv, 1, block_t, size - 1, comm);
}

/**
 * Загрузить входную конфигурацию.
 * Формат файла, число шагов, как часто сохранять, размер поля, затем идут
 * координаты заполненых клеток: steps save_steps nx ny i1 j2 i2 j2
 */
void life_init(const char *path, life_t *l) {
    FILE *fd = fopen(path, "r");
    assert(fd);
    assert(fscanf(fd, "%d\n", &l->steps));
    assert(fscanf(fd, "%d\n", &l->save_steps));
    // printf("Steps %d, save every %d step.\n", l->steps, l->save_steps);
    assert(fscanf(fd, "%d %d\n", &l->nx, &l->ny));
    // printf("Field size: %dx%d\n", l->nx, l->ny);

    l->u0 = (int *) calloc(l->nx * l->ny, sizeof(int));
    l->u1 = (int *) calloc(l->nx * l->ny, sizeof(int));

    int i, j, r, cnt;
    cnt = 0;
    while ((r = fscanf(fd, "%d %d\n", &i, &j)) != EOF) {
        l->u0[ind(i, j)] = 1;
        cnt++;
    }
    // printf("Loaded %d life cells.\n", cnt);
    fclose(fd);

    /* MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &(l->size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(l->rank));

    l->dims[0] = l->dims[1] = 0;
    MPI_Dims_create(l->size, 2, l->dims);
    unsigned int periods[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, l->dims, periods, 0, &(l->comm_cart));
    MPI_Cart_coords(l->comm_cart, l->rank, 2, l->coords);

    decomposition(l->nx, l->dims[0], l->coords[0], l->start, l->stop);
    decomposition(l->ny, l->dims[1], l->coords[1], l->start + 1, l->stop + 1);

    MPI_Comm_split(l->comm_cart, l->coords[1], l->coords[0], l->comm_dims);
    MPI_Comm_split(l->comm_cart, l->coords[0], l->coords[1], l->comm_dims + 1);

    int start[2], stop[2];
    decomposition(l->nx, l->dims[0], 0, start, stop);
    decomposition(l->ny, l->dims[1], l->coords[1], start + 1, stop + 1);

    MPI_Datatype block_raw;
    MPI_Type_vector(stop[1] - start[1], stop[0] - start[0], l->nx, MPI_INT, &block_raw);
    MPI_Type_create_resized(block_raw, 0, (stop[0] - start[0]) * sizeof(int), &(l->block_t));
    MPI_Type_commit(&(l->block_t));

    decomposition(l->ny, l->dims[1], 0, start + 1, stop + 1);
    MPI_Type_contiguous(l->nx * (stop[1] - start[1]), MPI_INT, &(l->row_t));
    MPI_Type_commit(&(l->row_t));

    MPI_Type_vector(1, l->stop[0] - l->start[0], 0, MPI_INT, &(l->row));
    MPI_Type_commit(&(l->row));

    MPI_Type_vector(l->stop[1] - l->start[1], 1, l->nx, MPI_INT, &(l->column));
    MPI_Type_commit(&(l->column));
}

void life_free(life_t *l) {
    free(l->u0);
    free(l->u1);
    l->nx = l->ny = 0;

    MPI_Type_free(&(l->block_t));
    MPI_Type_free(&(l->row_t));
    MPI_Type_free(&(l->row));
    MPI_Type_free(&(l->column));

    MPI_Comm_free(&(l->comm_cart));
    MPI_Comm_free(&(l->comm_dims[0]));
    MPI_Comm_free(&(l->comm_dims[1]));
}

void life_save_vtk(const char *path, life_t *l) {
    FILE *f;
    int i1, i2, j;

    struct stat st = {0};
    if (stat("vtk", &st) == -1) {
        mkdir("vtk", 0700);
    }

    f = fopen(path, "w");
    assert(f);
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Created by write_to_vtk2d\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET STRUCTURED_POINTS\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", l->nx + 1, l->ny + 1);
    fprintf(f, "SPACING %d %d 0.0\n", 1, 1);
    fprintf(f, "ORIGIN %d %d 0.0\n", 0, 0);
    fprintf(f, "CELL_DATA %d\n", l->nx * l->ny);

    fprintf(f, "SCALARS life int 1\n");
    fprintf(f, "LOOKUP_TABLE life_table\n");
    for (i2 = 0; i2 < l->ny; i2++) {
        for (i1 = 0; i1 < l->nx; i1++) {
            fprintf(f, "%d\n", l->u0[ind(i1, i2)]);
        }
    }
    fclose(f);
}

void life_step(life_t *l) {
    int i, j;
    for (j = l->start[1]; j < l->stop[1]; j++) {
        for (i = l->start[0]; i < l->stop[0]; i++) {
            int n = 0;
            n += l->u0[ind(i + 1, j)];
            n += l->u0[ind(i + 1, j + 1)];
            n += l->u0[ind(i, j + 1)];
            n += l->u0[ind(i - 1, j)];
            n += l->u0[ind(i - 1, j - 1)];
            n += l->u0[ind(i, j - 1)];
            n += l->u0[ind(i - 1, j + 1)];
            n += l->u0[ind(i + 1, j - 1)];
            l->u1[ind(i, j)] = 0;
            if (n == 3 && l->u0[ind(i, j)] == 0) {
                l->u1[ind(i, j)] = 1;
            }
            if ((n == 3 || n == 2) && l->u0[ind(i, j)] == 1) {
                l->u1[ind(i, j)] = 1;
            }
        }
    }
    int *tmp;
    tmp = l->u0;
    l->u0 = l->u1;
    l->u1 = tmp;
}

void decomposition(const int n, const int p, const int k, int *start, int *stop) {
    int l = n / p; // длинна куска
    *start = l * k;
    *stop = *start + l;
    if (k == p - 1)
        *stop = n;
}

void exchange_rows(life_t *l) {
    int this_stop = ind(l->start[0], l->stop[1] - 1);
    int left_stop = ind(l->start[0], l->start[1] - 1);
    int this_start = ind(l->start[0], l->start[1]);
    int right_start = ind(l->start[0], l->stop[1]);

    int right = (l->coords[1] + 1) % l->dims[1];
    int left = (l->coords[1] - 1 + l->dims[1]) % l->dims[1];

    MPI_Send(l->u0 + this_stop, 1, l->row, right, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + left_stop, 1, l->row, left, 0, l->comm_dims[1], MPI_STATUS_IGNORE);
    MPI_Send(l->u0 + this_start, 1, l->row, left, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + right_start, 1, l->row, right, 0, l->comm_dims[1], MPI_STATUS_IGNORE);
}

void exchange_columns(life_t *l) {
    int this_stop = ind(l->stop[0] - 1, l->start[1]);
    int left_stop = ind(l->start[0] - 1, l->start[1]);
    int this_start = ind(l->start[0], l->start[1]);
    int right_start = ind(l->stop[0], l->start[1]);

    int right = (l->coords[0] + 1) % l->dims[0];
    int left = (l->coords[0] - 1 + l->dims[0]) % l->dims[0];

    MPI_Send(l->u0 + this_stop, 1, l->column, right, 0, l->comm_dims[0]);
    MPI_Recv(l->u0 + left_stop, 1, l->column, left, 0, l->comm_dims[0], MPI_STATUS_IGNORE);
    MPI_Send(l->u0 + this_start, 1, l->column, left, 0, l->comm_dims[0]);
    MPI_Recv(l->u0 + right_start, 1, l->column, right, 0, l->comm_dims[0], MPI_STATUS_IGNORE);
}

void exchange_corners(life_t *l) {
    int right = (l->coords[1] + 1) % l->dims[1];
    int left = (l->coords[1] - 1 + l->dims[1]) % l->dims[1];

    MPI_Send(l->u0 + ind(l->stop[0], l->stop[1] - 1), 1, MPI_INT, right, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + ind(l->stop[0], l->start[1] - 1), 1, MPI_INT, left, 0, l->comm_dims[1], MPI_STATUS_IGNORE);

    MPI_Send(l->u0 + ind(l->start[0] - 1, l->stop[1] - 1), 1, MPI_INT, right, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + ind(l->start[0] - 1, l->start[1] - 1), 1, MPI_INT, left, 0, l->comm_dims[1], MPI_STATUS_IGNORE);

    MPI_Send(l->u0 + ind(l->start[0] - 1, l->start[1]), 1, MPI_INT, left, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + ind(l->start[0] - 1, l->stop[1]), 1, MPI_INT, right, 0, l->comm_dims[1], MPI_STATUS_IGNORE);

    MPI_Send(l->u0 + ind(l->stop[0], l->start[1]), 1, MPI_INT, left, 0, l->comm_dims[1]);
    MPI_Recv(l->u0 + ind(l->stop[0], l->stop[1]), 1, MPI_INT, right, 0, l->comm_dims[1], MPI_STATUS_IGNORE);
}

void life_exchange(life_t *l) {
    exchange_columns(l);
    exchange_rows(l);
    exchange_corners(l);
}