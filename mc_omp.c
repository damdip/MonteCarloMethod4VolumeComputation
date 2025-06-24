#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int is_inside(double *point, int dim); // Da implementare

double montecarlo_volume(int dim, int num_points, double *min, double *max) {
    int count_inside = 0;

    #pragma omp parallel
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num(); // Seed per random
        int local_count = 0;

        #pragma omp for
        for (int i = 0; i < num_points; i++) {
            double point[dim];
            for (int d = 0; d < dim; d++)
                point[d] = min[d] + ((double)rand_r(&seed) / RAND_MAX) * (max[d] - min[d]);

            if (is_inside(point, dim))
                local_count++;
        }

        #pragma omp atomic
        count_inside += local_count;
    }

    // Calcola volume del bounding box
    double volume_box = 1.0;
    for (int i = 0; i < dim; i++)
        volume_box *= (max[i] - min[i]);

    return (count_inside / (double)num_points) * volume_box;
}
