/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, 2016
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/stat.h>

#define ind(i, j) (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx))

typedef struct {
    int nx, ny;
    int *u0;
    int *u1;
    int steps;
    int save_steps;

    /* MPI */
    int start, stop;
    int rank;
    int size;
    MPI_Datatype block_type;
    MPI_Datatype string_type;
} life_t;

void life_init(const char *path, life_t *l);
void life_free(life_t *l);
void life_step(life_t *l);
void life_save_vtk(const char *path, life_t *l);
void life_collect(life_t *l);
void life_exchange(life_t *life_collect);
void decomposition(const int n, const int p, const int k, int *start, int *stop);

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

    life_free(&l);
    MPI_Finalize();
    return 0;
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

    l->u0 = (int *)calloc(l->nx * l->ny, sizeof(int));
    l->u1 = (int *)calloc(l->nx * l->ny, sizeof(int));

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
    decomposition(l->ny, l->size, l->rank, &(l->start), &(l->stop));

    int s1, s2;
    decomposition(l->ny, l->size, 0, &s1, &s2);
    // MPI_Type_vector(l->nx, s2 - s1, l->nx, MPI_INT, &(l->block_type));
    MPI_Type_contiguous(l->nx * (s2 - s1), MPI_INT, &(l->block_type));
    MPI_Type_contiguous(l->nx, MPI_INT, &(l->string_type));
    MPI_Type_commit(&(l->block_type));
    MPI_Type_commit(&(l->string_type));
}

void life_free(life_t *l) {
    free(l->u0);
    free(l->u1);
    l->nx = l->ny = 0;
    MPI_Type_free(&(l->block_type));
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
    for (j = l->start; j < l->stop; j++) {
        for (i = 0; i < l->nx; i++) {
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
    int l = n / p;  // длинна куска
    *start = l * k;
    *stop = *start + l;
    if (k == p - 1) *stop = n;
}

void life_collect(life_t *l) {
    if (l->rank == l->size - 1) {
        int i;
        for (i = 0; i < l->size - 1; i++) {
            int s1, s2;
            decomposition(l->ny, l->size, i, &s1, &s2);
            MPI_Recv(l->u0 + ind(0, s1), 1, l->block_type, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Send(l->u0 + ind(0, l->start), 1, l->block_type, l->size - 1, 0, MPI_COMM_WORLD);
    }
}

void life_exchange(life_t *l) {
    int left = (l->rank == 0) ? l->size - 1 : l->rank - 1;   // Левый сосед
    int right = (l->rank == l->size - 1) ? 0 : l->rank + 1;  // Правый сосед

    MPI_Send(l->u0 + ind(0, l->stop - 1), 1, l->string_type, left, 0, MPI_COMM_WORLD);
    MPI_Recv(l->u0 + ind(0, l->start - 1), 1, l->string_type, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Send(l->u0 + ind(0, l->start), 1, l->string_type, right, 0, MPI_COMM_WORLD);
    MPI_Recv(l->u0 + ind(0, l->stop), 1, l->string_type, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}