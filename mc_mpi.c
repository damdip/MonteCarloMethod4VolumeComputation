#include <mpi.h>

// Ogni processo genera parte dei punti e somma i risultati
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Parametri globali
    int dim = 3, num_points = 1000000;
    double min[3] = {0,0,0}, max[3] = {1,1,1}; // esempio 3D

    int local_N = num_points / size;
    int local_inside = 0;

    // Simulazione
    srand(time(NULL) + rank);
    for (int i = 0; i < local_N; i++) {
        double point[3];
        for (int d = 0; d < dim; d++)
            point[d] = min[d] + ((double)rand() / RAND_MAX) * (max[d] - min[d]);

        if (is_inside(point, dim))
            local_inside++;
    }

    // Riduzione: somma dei risultati locali
    int total_inside = 0;
    MPI_Reduce(&local_inside, &total_inside, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double volume_box = 1.0;
        for (int i = 0; i < dim; i++) volume_box *= (max[i] - min[i]);
        double volume = (total_inside / (double)num_points) * volume_box;
        printf("Stima volume: %f\n", volume);
    }

    MPI_Finalize();
    return 0;
}
