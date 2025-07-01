#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Struttura per definire i limiti di ogni dimensione
typedef struct {
    double min;
    double max;
} Bounds;

// Struttura per un punto k-dimensionale
typedef struct {
    double *coords;
    int dimension;
} Point;

// Funzione per creare un punto k-dimensionale
Point* create_point(int k) {
    Point *p = malloc(sizeof(Point));
    p->coords = malloc(k * sizeof(double));
    p->dimension = k;
    return p;
}

// Funzione per liberare la memoria di un punto
void free_point(Point *p) {
    free(p->coords);
    free(p);
}

// Genera un numero random tra min e max
double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Genera un punto casuale nell'ipercubo k-dimensionale
Point* generate_random_point(Bounds *bounds, int k) {
    Point *p = create_point(k);
    for (int i = 0; i < k; i++) {
        p->coords[i] = random_double(bounds[i].min, bounds[i].max);
    }
    return p;
}

// Calcola il volume dell'ipercubo k-dimensionale
double calculate_hypercube_volume(Bounds *bounds, int k) {
    double volume = 1.0;
    for (int i = 0; i < k; i++) {
        volume *= (bounds[i].max - bounds[i].min);
    }
    return volume;
}


// Ipersfera k-dimensionale centrata nell'origine
int is_inside_hypersphere(Point *p, double radius) {
    double sum_squares = 0.0;
    for (int i = 0; i < p->dimension; i++) {
        sum_squares += p->coords[i] * p->coords[i];
    }
    return sum_squares <= radius * radius;
}


// Funzione principale per il metodo Monte Carlo
double monte_carlo_volume(Bounds *bounds, int k, long n_samples, 
                         int (*shape_function)(Point*)) {
    
    srand(time(NULL)); // Inizializza il generatore di numeri casuali
    
    long points_inside = 0;
    double hypercube_volume = calculate_hypercube_volume(bounds, k);
    
    printf("Inizio simulazione Monte Carlo...\n");
    printf("Dimensioni: %d\n", k);
    printf("Campioni: %ld\n", n_samples);
    printf("Volume ipercubo: %.6f\n\n", hypercube_volume);
    
    // Genera punti e conta quelli dentro la figura
    for (long i = 0; i < n_samples; i++) {
        Point *p = generate_random_point(bounds, k);
        
        if (shape_function(p)) {
            points_inside++;
        }
        
        free_point(p);
        
        /* Progresso ogni 10% dei campioni
        if (i % (n_samples / 10) == 0 && i > 0) {
            double current_estimate = hypercube_volume * 
                                    ((double)points_inside / i);
            printf("Progresso: %.1f%% - Stima corrente: %.6f\n", 
                   (double)i / n_samples * 100, current_estimate);
        }
        */

    }

    //in versione open mp prendere l'array dove sono stati messi i "points inside" o fare reduction
    // in versione open mpi fare una reduction sicuramente. 
    
    double estimated_volume = hypercube_volume * 
                             ((double)points_inside / n_samples);
    
    printf("\nRisultati finali:\n");
    printf("Punti totali: %ld\n", n_samples);
    printf("Punti dentro la figura: %ld\n", points_inside);
    printf("Rapporto: %.6f\n", (double)points_inside / n_samples);
    printf("Volume stimato: %.6f\n", estimated_volume);
    
    return estimated_volume;
}

// Esempio di uso per ipersfera
double hypersphere_wrapper(Point *p) {
    return is_inside_hypersphere(p, 1.0); // Raggio = 1
}

int main(int argc , char** argv) {
    //Si passa da linea di comando il parametro k
    int k;
    if( argc > 1){
        k = atoi(argv[1]);
    }else{
        k = 3;    
    }
    printf("=== METODO MONTE CARLO PER VOLUMI %d-DIMENSIONALI ===\n\n", k);
    
    // ESEMPIO 1: Ipersfera 3D con raggio 1
    printf("ESEMPIO: Ipersfera %dD (raggio = 1)\n", k);
    float volume_effettivo = pow(M_PI , k/2.0) / tgamma(k/2.0 + 1.0);
    printf("Volume teorico: %.6f\n", volume_effettivo); // 4/3 * π
    


    // Creazione bounds
    Bounds bounds[k];

    for ( int i = 0 ; i< k ; i++ ){
        bounds[i].min = -2;
        bounds[i].max = +2;
        
    }
    
    double volume_sphere = monte_carlo_volume(bounds, k, 100000000, 
                                            hypersphere_wrapper);
    
    printf("Errore relativo: %.2f%%\n", 
           fabs(volume_sphere - volume_effettivo) / (volume_effettivo) * 100);
    
    return 0;
}

